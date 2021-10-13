
rule mutect2_single:
    input: tumor=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam")
    output: vcf=os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{samples}.{chroms}.vcf"),
            # vcfgz=os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{samples}.{chroms}.vcf.gz"),
            read_orientation_file=os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{samples}.{chroms}.f1r2.tar.gz"),
            statsfiles=os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{samples}.{chroms}.vcf.stats")
    params: tumorsample="{samples}",
            genome = config['references']['GENOME'],
            pon=config['references']['PON'],
            germsource=config['references']['GNOMAD'],
            ver_gatk=config['tools']['gatk4']['version'],
            # ver_bcftools=config['tools']['bcftools']['version'],
            rname="mutect2"
    threads: 2
    shell: """
           if [ ! -d "$(dirname {output.vcf})" ]; then
               mkdir -p "$(dirname {output.vcf})"
           fi
           module load GATK/{params.ver_gatk}
           gatk Mutect2 -R {params.genome} -I {input.tumor} --panel-of-normals {params.pon} --germline-resource {params.germsource} -L {wildcards.chroms} -O {output.vcf} --f1r2-tar-gz {output.read_orientation_file} --independent-mates
           """

        
rule pileup_single:
    input: tumor=os.path.join(output_bamdir,"final_bams","{samples}.bam"),
           bai=os.path.join(output_bamdir,"final_bams","{samples}.bai"),
           intervals=intervals_file
    output: pileup=temp(os.path.join(output_somatic_snpindels,"mutect2_out","pileup_summaries","{samples}.pileup.table")),
    params: genome=config['references']['GENOME'],germsource=config['references']['1000GSNP'],ver_gatk=config['tools']['gatk4']['version'],chroms=chroms, rname="pileup"
    shell: """
           module load GATK/{params.ver_gatk}
           gatk --java-options "-Xmx10g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID" GetPileupSummaries -R {params.genome} -I {input.tumor} -V {params.germsource} -L {input.intervals} -O {output.pileup}
           """
           
rule contamination_single:
    input: pileup=os.path.join(output_somatic_snpindels,"mutect2_out","pileup_summaries","{samples}.pileup.table")
    output: tumor_summary=os.path.join(output_somatic_base,"qc","gatk_contamination","{samples}.contamination.table")
    # resources: tmpdir="/lscratch/$SLURM_JOB_ID"
    params: genome=config['references']['GENOME'],germsource=config['references']['1000GSNP'],ver_gatk=config['tools']['gatk4']['version'],chroms=chroms, rname="pileup"
    # group: "contam"
    shell: """
           module load GATK/{params.ver_gatk}
           gatk CalculateContamination -I {input.pileup} -O {output.tumor_summary}
           """

rule mutect_single:
       input:  tumor=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam"),
       output: vcf=os.path.join(output_somatic_snpindels,"mutect_out","chrom_split","{samples}.{chroms}.vcf"),
               # vcfgz=temp(os.path.join(output_somatic_snpindels,"mutect_out","chrom_split","{samples}.{chroms}.vcf.gz")),
               stats=os.path.join(output_somatic_snpindels,"mutect_out","chrom_split","{samples}.{chroms}.stats.out"),
               # final=os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.FINAL.vcf"),
       params: genome=config['references']['GENOME'],
               pon=config['references']['PON'],
               cosmic=config['references']['COSMIC'],
               dbsnp=config['references']['DBSNP'],
               ver_mutect=config['tools']['mutect']['version'],
               # ver_bcftools=config['tools']['bcftools']['version'],
               rname="mutect"
       shell:  """
               myoutdir="$(dirname {output.vcf})"
               if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
               
               module load muTect/{params.ver_mutect}
               muTect --analysis_type MuTect --reference_sequence {params.genome} --normal_panel {params.pon} --vcf {output.vcf} --cosmic {params.cosmic} --dbsnp {params.dbsnp} -L {wildcards.chroms} --disable_auto_index_creation_and_locking_when_reading_rods --input_file:tumor {input.tumor} --out {output.stats} -rf BadCigar
               """

rule mutect_filter_single:
    input: vcf=os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.collected.vcf"),
    output: final=os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.FINAL.vcf"),
            norm=os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.FINAL.norm.vcf"),
    params: tumorsample="{samples}",
            genome=config['references']['GENOME'],
            ver_gatk=config['tools']['gatk4']['version'],
            ver_bcftools=config['tools']['bcftools']['version'],
            rname="mutect_filter"
    shell: """
           module load GATK/{params.ver_gatk}
           gatk SelectVariants -R {params.genome} --variant {input.vcf} --exclude-filtered --output {output.final}
           
           module load bcftools/{params.ver_bcftools}
           ## (At least) VarScan outputs ambiguous IUPAC bases/codes sometimes; the piped awk one-liner sets them to N
           ## From: https://github.com/fpbarthel/GLASS/issues/23
           bcftools sort -T /lscratch/$SLURM_JOB_ID "{output.final}" | \
               bcftools norm --threads $((SLURM_CPUS_PER_TASK - 1)) --check-ref s -f {params.genome} -O v | \
               awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' | sed '/^$/d' > {output.norm}
           """

rule vardict_single:
       input:  tumor=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam"),
       output: vcf=os.path.join(output_somatic_snpindels,"vardict_out","chrom_split","{samples}.{chroms}.vcf"),
               # vcfgz=temp(os.path.join(output_somatic_snpindels,"vardict_out","chrom_split","{samples}.{chroms}.vcf.gz")),
       params: genome=config['references']['GENOME'],
               targets=exome_targets_bed,
               pon=config['references']['PON'],
               ver_bcftools=config['tools']['bcftools']['version'],
               rname="vardict"
       shell:  """
               myoutdir="$(dirname {output.vcf})"
               if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
               module load R/3.6.1
               module load samtools/1.8
               #module load java/1.8.0_92
               
               /data/CCBR_Pipeliner/db/PipeDB/bin/VarDictJava/build/install/VarDict/bin/VarDict -G {params.genome} -f 0.05 -x 500 --nosv -b {input.tumor} -t -Q 20 -c 1 -S 2 -E 3 {params.targets} | \
               /data/CCBR_Pipeliner/db/PipeDB/bin/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
               /data/CCBR_Pipeliner/db/PipeDB/bin/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl -N {wildcards.samples} -Q 20 -d 10 -v 6 -S -E -f 0.05 > {output.vcf}
               """

rule vardict_filter_single:
       input:  vcf=os.path.join(output_somatic_snpindels,"vardict_out","vcf","{samples}.collected.vcf"),
       output: final=os.path.join(output_somatic_snpindels,"vardict_out","vcf","{samples}.FINAL.vcf"),
               norm=os.path.join(output_somatic_snpindels,"vardict_out","vcf","{samples}.FINAL.norm.vcf"),
               # filtered=temp(os.path.join(output_somatic_snpindels,"vardict_out","vcf","{samples}.filtered.vcf")),
       params: tumorsample="{samples}",
               genome=config['references']['GENOME'],
               targets=exome_targets_bed,pon=config['references']['PON'],
               ver_gatk=config['tools']['gatk4']['version'],
               ver_bcftools=config['tools']['bcftools']['version'],
               rname="vardict_filter"
       shell:  """
               module load GATK/{params.ver_gatk}
               gatk SelectVariants -R {params.genome} --variant {input.vcf} --discordance {params.pon} --exclude-filtered --output {output.final}
               
               module load bcftools/{params.ver_bcftools}
               ## (At least) VarScan outputs ambiguous IUPAC bases/codes sometimes; the piped awk one-liner sets them to N
               ## From: https://github.com/fpbarthel/GLASS/issues/23
               bcftools sort -T /lscratch/$SLURM_JOB_ID "{output.final}" | \
                   bcftools norm --threads $((SLURM_CPUS_PER_TASK - 1)) --check-ref s -f {params.genome} -O v | \
                   awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' | sed '/^$/d' > {output.norm}
               """

rule varscan_single:
       input:  tumor=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam"),
       output: vcf=os.path.join(output_somatic_snpindels,"varscan_out","chrom_split","{samples}.{chroms}.vcf"),
       params: genome=config['references']['GENOME'],
               ver_varscan=config['tools']['varscan']['version'],
               ver_bcftools=config['tools']['bcftools']['version'],
               rname="varscan"
       shell:  """
               if [ ! -d "$(dirname {output.vcf})" ]; then
                   mkdir -p "$(dirname {output.vcf})"
               fi
               module load VarScan/{params.ver_varscan}
               
               varscan_opts="--strand-filter 0 --min-var-freq 0.01 --output-vcf 1 --variants 1"
               
               pileup_cmd="samtools mpileup -d 100000 -q 15 -Q 15 -f {params.genome} {input.tumor}"
               varscan_cmd="varscan mpileup2cns <($pileup_cmd) $varscan_opts"
               
               #full_cmd="$varscan_cmd | bcftools view -U > {output.vcf}"
               #eval "$full_cmd"
               
               eval "$varscan_cmd > {output.vcf}.gz"
               eval "bcftools view -U {output.vcf}.gz > {output.vcf}"
               
               """

rule varscan_filter_single:
    input: vcf=os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.collected.vcf"),
    output: filtered=temp(os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.filtered.vcf")),
            filtered1=temp(os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.filtered.1.vcf")),
            samplesfile=temp(os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.FINAL.vcf.samples")),
            final=os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.FINAL.vcf"),
            norm=os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.FINAL.norm.vcf"),
    params: tumorsample="{samples}",genome=config['references']['GENOME'],
            pon=config['references']['PON'],
            basedir=BASEDIR,
            ver_gatk=config['tools']['gatk4']['version'],
            ver_varscan=config['tools']['varscan']['version'],
            ver_bcftools=config['tools']['bcftools']['version'],
            rname="varscan_filter"
    shell: """
           module load VarScan/{params.ver_varscan}
           varscan filter {input.vcf} {config[tools][varscan][filter_settings]} > {output.filtered1}
           
           module load GATK/{params.ver_gatk}
           gatk SelectVariants -R {params.genome} --variant {output.filtered1} --discordance {params.pon} --exclude-filtered --output {output.filtered}
           samplesFile="{output.samplesfile}"
           echo -e "TUMOR\t{params.tumorsample}\n" > "{output.samplesfile}"
           
           module load bcftools/{params.ver_bcftools}
           bcftools reheader -o "{output.final}" -s "{output.samplesfile}" "{output.filtered}"
           ## (At least) VarScan outputs ambiguous IUPAC bases/codes sometimes; the piped awk one-liner sets them to N
           ## From: https://github.com/fpbarthel/GLASS/issues/23
           bcftools sort -T /lscratch/$SLURM_JOB_ID "{output.final}" | \
            bcftools norm --threads $((SLURM_CPUS_PER_TASK - 1)) --check-ref s -f {params.genome} -O v | \
            awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' | sed '/^$/d' > {output.norm}
           """
           

