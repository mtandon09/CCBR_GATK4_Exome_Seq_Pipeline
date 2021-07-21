
rule split_bam_by_chrom:
      input:  bam=os.path.join(output_bamdir,"final_bams","{samples}.bam"),
              bai=os.path.join(output_bamdir,"final_bams","{samples}.bam.bai"),
      output: split_bam=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam"),
              split_bam_idx=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bai")
      params: ver_samtools=config['tools']['samtools']['version'],rname="bam_split"
      shell:  """
              if [ ! -d "$(dirname {output.split_bam})" ]; then
                mkdir -p "$(dirname {output.split_bam})"
              fi
              module load samtools/{params.ver_samtools}
              samtools view -b -o {output.split_bam} -@ $((SLURM_CPUS_PER_TASK-1)) {input.bam} {wildcards.chroms}
              samtools index -@ 2 {output.split_bam} {output.split_bam_idx}
              cp {output.split_bam_idx} {output.split_bam}.bai
              """

rule gatk_mutect2:
    input: normal = lambda w: [os.path.join(output_bamdir,"chrom_split",pairs_dict[w.samples] + ".{chroms}.split.bam")],
           tumor=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam")
    output: vcf=os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{samples}.{chroms}.vcf"),
            read_orientation_file=os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{samples}.{chroms}.f1r2.tar.gz"),
            statsfiles=os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{samples}.{chroms}.vcf.stats")
    params: normalsample=lambda w: [pairs_dict[w.samples]],
            tumorsample="{samples}",
            chrom="{chroms}",
            genome = config['references']['GENOME'],
            pon=config['references']['PON'],
            germsource=config['references']['GNOMAD'],
            ver_gatk=config['tools']['gatk4']['version'],
            rname="mutect2"
    threads: 2
    shell: """
           if [ ! -d "$(dirname {output.vcf})" ]; then
               mkdir -p "$(dirname {output.vcf})"
           fi
           module load GATK/{params.ver_gatk}
           gatk Mutect2 -R {params.genome} -I {input.tumor} -I {input.normal} -normal {params.normalsample} --panel-of-normals {params.pon} --germline-resource {params.germsource} -L {params.chrom} -O {output.vcf} --f1r2-tar-gz {output.read_orientation_file} --independent-mates
           """

rule LearnReadOrientationModel:
    input: vcf=expand(os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{{samples}}.{chroms}.vcf"), chroms=chroms),
           read_orientation_file=expand(os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{{samples}}.{chroms}.f1r2.tar.gz"), chroms=chroms)
    output: model=os.path.join(output_somatic_snpindels,"mutect2_out","read_orientation_data","{samples}.read-orientation-model.tar.gz")
    params: normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",genome = config['references']['GENOME'],ver_gatk=config['tools']['gatk4']['version'],rname="LearnReadOrientationModel"
    shell: """
           module load GATK/{params.ver_gatk}
           input_str="--input $(echo "{input.read_orientation_file}" | sed -e 's/ / --input /g')"
           
           gatk LearnReadOrientationModel --output {output.model} $input_str
           """

rule pileup_paired:
    input: tumor=os.path.join(output_bamdir,"final_bams","{samples}.bam"),
           tumorbai=os.path.join(output_bamdir,"final_bams","{samples}.bai"),
           normal = lambda w: [os.path.join(output_bamdir,"final_bams", pairs_dict[w.samples] + ".bam")],
           normalbai=lambda w: [os.path.join(output_bamdir,"final_bams", pairs_dict[w.samples] + ".bai")],
           intervals=intervals_file
    output: tumor_summary=os.path.join(output_somatic_snpindels,"mutect2_out","pileup_summaries","{samples}_tumor.pileup.table"),
            normal_summary=os.path.join(output_somatic_snpindels,"mutect2_out","pileup_summaries","{samples}_normal.pileup.table")
    params: genome=config['references']['GENOME'],germsource=config['references']['1000GSNP'],ver_gatk=config['tools']['gatk4']['version'],rname="pileup"
    shell: """
           module load GATK/{params.ver_gatk}
           gatk --java-options '-Xmx48g' GetPileupSummaries -I {input.tumor} -V {params.germsource} -L {input.intervals} -O {output.tumor_summary} & \
           gatk --java-options '-Xmx48g' GetPileupSummaries -I {input.normal} -V {params.germsource} -L {input.intervals} -O {output.normal_summary} & \
           wait
           """

# rule pileup_normal:
#     input: normal = lambda w: [os.path.join(output_bamdir, pairs_dict[w.samples] + ".recal.bam")],
#     output: summary=os.path.join(output_somatic_snpindels,"mutect2_out","pileup_summaries","{samples}_normal.pileup.table")
#     params: genome=config['references']['GENOME'],germsource=config['references']['1000GSNP'],ver_gatk=config['tool_versions']['gatk4'],rname="pileup"
#     shell: """
#            module load GATK/{params.ver_gatk}
#            gatk --java-options '-Xmx48g' GetPileupSummaries -I {input.normal} -V {params.germsource} -L {params.germsource} -O {output.summary}
#            """

rule contamination_paired:
    input: tumor=os.path.join(output_somatic_snpindels,"mutect2_out","pileup_summaries","{samples}_tumor.pileup.table"),
           normal=os.path.join(output_somatic_snpindels,"mutect2_out","pileup_summaries","{samples}_normal.pileup.table"),
    output: tumor_summary=os.path.join(output_somatic_base,"qc","gatk_contamination","{samples}.contamination.table"),
    params: normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",genome=config['references']['GENOME'],ver_gatk=config['tools']['gatk4']['version'],rname="contamination"
    shell: """
           module load GATK/{params.ver_gatk}
           gatk CalculateContamination -I {input.tumor} --matched-normal {input.normal} -O {output.tumor_summary}
           """

rule mutect2_filter:
    input: vcf=os.path.join(output_somatic_snpindels,"mutect2_out","vcf","{samples}.collected.vcf"),
           summary=os.path.join(output_somatic_base,"qc","gatk_contamination","{samples}.contamination.table"),
           model=os.path.join(output_somatic_snpindels,"mutect2_out","read_orientation_data","{samples}.read-orientation-model.tar.gz"),
           statsfiles=expand(os.path.join(output_somatic_snpindels,"mutect2_out","chrom_split","{{samples}}.{chroms}.vcf.stats"), chroms=chroms)
    output: marked_vcf=os.path.join(output_somatic_snpindels,"mutect2_out","vcf","{samples}.filtered.vcf"),
            final=os.path.join(output_somatic_snpindels,"mutect2_out","vcf","{samples}.FINAL.vcf"),
            norm=os.path.join(output_somatic_snpindels,"mutect2_out","vcf","{samples}.FINAL.norm.vcf"),
    params: normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",genome=config['references']['GENOME'],
            ver_gatk=config['tools']['gatk4']['version'],
            ver_bcftools=config['tools']['bcftools']['version'],
            rname="mutect2_filter"
    shell: """
           module load GATK/{params.ver_gatk}
           
           statfiles="--stats $(echo "{input.statsfiles}" | sed -e 's/ / --stats /g')"
           gatk MergeMutectStats $statfiles -O {output.final}.stats
           gatk FilterMutectCalls -R {params.genome} -V {input.vcf} --ob-priors {input.model} --contamination-table {input.summary} -O {output.marked_vcf} --stats {output.final}.stats
           gatk SelectVariants -R {params.genome} --variant {output.marked_vcf} --exclude-filtered --output {output.final}
           
           module load bcftools/{params.ver_bcftools}
           ## (At least) VarScan outputs ambiguous IUPAC bases/codes sometimes; the piped awk one-liner sets them to N
           ## From: https://github.com/fpbarthel/GLASS/issues/23
           bcftools sort -T /lscratch/$SLURM_JOB_ID "{output.final}" | \
               bcftools norm --threads $((SLURM_CPUS_PER_TASK - 1)) --check-ref s -f {params.genome} -O v | \
               awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' | sed '/^$/d' > {output.norm}
           """

           
rule strelka:
    input: normal = lambda w: [os.path.join(output_bamdir,"chrom_split",pairs_dict[w.samples] + ".{chroms}.split.bam")],
           tumor=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam")
    output: vcf=os.path.join(output_somatic_snpindels,"strelka_out","chrom_split","{samples}.{chroms}.vcf"),
            # final=os.path.join(output_somatic_snpindels,"strelka_out","chrom_split","{samples}.{chroms}.filtered.vcf"),
    params: genome=config['references']['GENOME'],
            pon=config['references']['PON'],
            basedir=BASEDIR,
            ver_strelka=config['tools']['strelka']['version'],
            rname="strelka"
    shell: """
           module load strelka/{params.ver_strelka}
           workdir={params.basedir}
           myoutdir="$(dirname {output.vcf})/{wildcards.samples}/{wildcards.chroms}"
           if [ -d "$myoutdir" ]; then rm -r "$myoutdir"; fi
           mkdir -p "$myoutdir"
           
           configureStrelkaSomaticWorkflow.py --ref={params.genome} --tumor={input.tumor} --normal={input.normal} --runDir="$myoutdir" --exome
           cd "$myoutdir"
           ./runWorkflow.py -m local -j $((SLURM_CPUS_PER_TASK-1))
           
           ##### I THINK STRELKA IS LOADING AN OLD VERSION OF JAVA CUZ GATK FAILS
           ##### 'module purge' fixes that on Biowulf, but for cross-architecture support, this should prob be split into a new rule
           
           module purge
           module load GATK/3.8-1
           GATK -m 12G CombineVariants -R {params.genome} --variant results/variants/somatic.snvs.vcf.gz --variant results/variants/somatic.indels.vcf.gz --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o "$(basename {output.vcf})"
           
           cd $workdir
           mv "$myoutdir/$(basename {output.vcf})" "{output.vcf}"
           
           
           """

rule strelka_filter:
    input: vcf=os.path.join(output_somatic_snpindels,"strelka_out","vcf","{samples}.collected.vcf"),
    output: filtered=temp(os.path.join(output_somatic_snpindels,"strelka_out","vcf","{samples}.filtered.vcf")),
            samplesfile=temp(os.path.join(output_somatic_snpindels,"strelka_out","vcf","{samples}.FINAL.vcf.samples")),
            final=os.path.join(output_somatic_snpindels,"strelka_out","vcf","{samples}.FINAL.vcf"),
            norm=os.path.join(output_somatic_snpindels,"strelka_out","vcf","{samples}.FINAL.norm.vcf"),
    params: normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",
            genome=config['references']['GENOME'],
            pon=config['references']['PON'],
            basedir=BASEDIR,
            ver_gatk=config['tools']['gatk4']['version'],
            ver_bcftools=config['tools']['bcftools']['version'],
            rname="strelka_filter"
    shell: """
           module load GATK/{params.ver_gatk}
           gatk SelectVariants -R {params.genome} --variant {input.vcf} --discordance {params.pon} --exclude-filtered --output {output.filtered}
           samplesFile="{output.samplesfile}"
           echo -e "TUMOR\t{params.tumorsample}\nNORMAL\t{params.normalsample}" > "{output.samplesfile}"
           
           echo "Reheading VCFs with sample names..."
           module load bcftools/{params.ver_bcftools}
           bcftools reheader -o "{output.final}" -s "{output.samplesfile}" "{output.filtered}"
           
           ## (At least) VarScan outputs ambiguous IUPAC bases/codes sometimes; the piped awk one-liner sets them to N
           ## From: https://github.com/fpbarthel/GLASS/issues/23
           bcftools sort -T /lscratch/$SLURM_JOB_ID "{output.final}" | \
               bcftools norm --threads $((SLURM_CPUS_PER_TASK - 1)) --check-ref s -f {params.genome} -O v | \
               awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' | sed '/^$/d' > {output.norm}
           """

rule mutect:
       input:  normal = lambda w: [os.path.join(output_bamdir,"chrom_split",pairs_dict[w.samples] + ".{chroms}.split.bam")],
               tumor=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam"),
       output: vcf=os.path.join(output_somatic_snpindels,"mutect_out","chrom_split","{samples}.{chroms}.vcf"),
               stats=os.path.join(output_somatic_snpindels,"mutect_out","chrom_split","{samples}.{chroms}.stats.out"),
               # final=os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.FINAL.vcf"),
       params: normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",pon=config['references']['PON'],genome=config['references']['GENOME'],cosmic=config['references']['COSMIC'],dbsnp=config['references']['DBSNP'],rname="mutect"
       shell:  """
               myoutdir="$(dirname {output.vcf})"
               if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
               
               module load muTect/1.1.7
               muTect --analysis_type MuTect --reference_sequence {params.genome} --normal_panel {params.pon} --vcf {output.vcf} --cosmic {params.cosmic} --dbsnp {params.dbsnp} --disable_auto_index_creation_and_locking_when_reading_rods --input_file:normal {input.normal} --input_file:tumor {input.tumor} --out {output.stats} -rf BadCigar
               """

rule mutect_filter:
    input: vcf=os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.collected.vcf"),
    output: final=os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.FINAL.vcf"),
            norm=os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.FINAL.norm.vcf"),
            # filtered=temp(os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.filtered.vcf")),
            # samplesfile=temp(os.path.join(output_somatic_snpindels,"mutect_out","vcf","{samples}.FINAL.vcf.samples")),
    params: normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",
            genome=config['references']['GENOME'],
            pon=config['references']['PON'],
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
           
rule vardict:
       input:  normal = lambda w: [os.path.join(output_bamdir,"chrom_split",pairs_dict[w.samples] + ".{chroms}.split.bam")],
               tumor=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam"),
       output: vcf=os.path.join(output_somatic_snpindels,"vardict_out","chrom_split","{samples}.{chroms}.vcf"),
       params: normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",genome=config['references']['GENOME'],targets=exome_targets_bed,pon=config['references']['PON'],rname="vardict"
       shell:  """
               myoutdir="$(dirname {output.vcf})"
               if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
               module load R/3.6.1
               module load samtools/1.8
               #module load java/1.8.0_92
               /data/CCBR_Pipeliner/db/PipeDB/bin/VarDictJava/build/install/VarDict/bin/VarDict -G {params.genome} -f 0.05 -N \"{params.tumorsample}|{params.normalsample}\" --nosv -b \"{input.tumor}|{input.normal}\" -t -Q 20 -c 1 -S 2 -E 3 {params.targets} | /data/CCBR_Pipeliner/db/PipeDB/bin/VarDictJava/build/install/VarDict/bin/testsomatic.R | /data/CCBR_Pipeliner/db/PipeDB/bin/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl -S -Q 20 -d 10 -M -N \"{params.tumorsample}|{params.normalsample}\" -f 0.05 > {output.vcf}
               """
rule vardict_filter:
       input:  vcf=os.path.join(output_somatic_snpindels,"vardict_out","vcf","{samples}.collected.vcf"),
       output: final=os.path.join(output_somatic_snpindels,"vardict_out","vcf","{samples}.FINAL.vcf"),
               filtered=temp(os.path.join(output_somatic_snpindels,"vardict_out","vcf","{samples}.filtered.vcf")),
               norm=os.path.join(output_somatic_snpindels,"vardict_out","vcf","{samples}.FINAL.norm.vcf"),
       params: normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",
               genome=config['references']['GENOME'],
               targets=exome_targets_bed,
               pon=config['references']['PON'],
               ver_gatk=config['tools']['gatk4']['version'],
               ver_bcftools=config['tools']['bcftools']['version'],
               rname="vardict"
       shell:  """
               module load bcftools
               bcftools filter --exclude \'STATUS=\"Germline\" | STATUS=\"LikelyLOH\" | STATUS=\"AFDiff\"\' {input.vcf} > {output.filtered}
               module load GATK/{params.ver_gatk}
               gatk SelectVariants -R {params.genome} --variant {output.filtered} --discordance {params.pon} --exclude-filtered --output {output.final}
               
               module load bcftools/{params.ver_bcftools}
               ## (At least) VarScan outputs ambiguous IUPAC bases/codes sometimes; the piped awk one-liner sets them to N
               ## From: https://github.com/fpbarthel/GLASS/issues/23
               bcftools sort -T /lscratch/$SLURM_JOB_ID "{output.final}" | \
               bcftools norm --threads $((SLURM_CPUS_PER_TASK - 1)) --check-ref s -f {params.genome} -O v | \
               awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' | sed '/^$/d' > {output.norm}
               """

rule varscan:
       input:  normal = lambda w: [os.path.join(output_bamdir,"chrom_split",pairs_dict[w.samples] + ".{chroms}.split.bam")],
               tumor=os.path.join(output_bamdir,"chrom_split","{samples}.{chroms}.split.bam"),
               tumor_summary=os.path.join(output_somatic_base,"qc","gatk_contamination","{samples}.tumor.contamination.table"),
               normal_summary=lambda w: [os.path.join(output_somatic_base,"qc","gatk_contamination","{samples}.normal.contamination.table")],
               # sequenza_purity=os.path.join(output_somatic_cnv,"sequenza_out","{samples}_alternative_solutions.txt")
       output: vcf=os.path.join(output_somatic_snpindels,"varscan_out","chrom_split","{samples}.{chroms}.vcf"),
               # samplesfile=os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.samples")
       params: genome=config['references']['GENOME'],normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",ver_varscan=config['tools']['varscan']['version'],rname="varscan"
       shell:  """
               if [ ! -d "$(dirname {output.vcf})" ]; then
                   mkdir -p "$(dirname {output.vcf})"
               fi
               module load VarScan/{params.ver_varscan}
               
               tumor_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 {input.tumor_summary} | cut -f2 ))" | bc -l)
               normal_purity=$( echo "1-$(printf '%.6f' $(tail -n -1 {input.normal_summary} | cut -f2 ))" | bc -l)
               
               varscan_opts="--strand-filter 1 --min-var-freq 0.01 --min-avg-qual 30 --somatic-p-value 0.05 --output-vcf 1 --normal-purity $normal_purity --tumor-purity $tumor_purity"
               
               dual_pileup="samtools mpileup -d 10000 -q 15 -Q 15 -f {params.genome} {input.normal} {input.tumor}"
               varscan_cmd="varscan somatic <($dual_pileup) {output.vcf} $varscan_opts --mpileup 1"
               
               eval "$varscan_cmd"
                      
               module load GATK/3.8-1
               GATK -m 12G CombineVariants -R {params.genome} --variant {output.vcf}.snp --variant {output.vcf}.indel --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {output.vcf}     
               """

rule varscan_filter:
    input: vcf=os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.collected.vcf"),
    output: filtered=temp(os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.filtered.vcf")),
            filtered1=temp(os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.filtered.1.vcf")),
            samplesfile=temp(os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.FINAL.vcf.samples")),
            final=os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.FINAL.vcf"),
            norm=os.path.join(output_somatic_snpindels,"varscan_out","vcf","{samples}.FINAL.norm.vcf"),
    params: normalsample=lambda w: [pairs_dict[w.samples]],tumorsample="{samples}",
            genome=config['references']['GENOME'],
            pon=config['references']['PON'],
            basedir=BASEDIR,
            filter_settings=config["tools"]["varscan"]["filter_settings"],
            ver_varscan=config['tools']['varscan']['version'],
            ver_gatk=config['tools']['gatk4']['version'],
            ver_bcftools=config['tools']['bcftools']['version'],
            rname="varscan_filter"
    shell: """
           module load VarScan/{params.ver_varscan}
           varscan filter {input.vcf} {params.filter_settings} > {output.filtered1}
           
           module load GATK/{params.ver_gatk}
           gatk SelectVariants -R {params.genome} --variant {output.filtered1} --discordance {params.pon} --exclude-filtered --output {output.filtered}
           samplesFile="{output.samplesfile}"
           echo -e "TUMOR\t{params.tumorsample}\nNORMAL\t{params.normalsample}" > "{output.samplesfile}"
           
           module load bcftools/{params.ver_bcftools}
           bcftools reheader -o "{output.final}" -s "{output.samplesfile}" "{output.filtered}"
           ## (At least) VarScan outputs ambiguous IUPAC bases/codes sometimes; the piped awk one-liner sets them to N
           ## From: https://github.com/fpbarthel/GLASS/issues/23
           bcftools sort -T /lscratch/$SLURM_JOB_ID "{output.final}" | \
            bcftools norm --threads $((SLURM_CPUS_PER_TASK - 1)) --check-ref s -f {params.genome} -O v | \
            awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' | sed '/^$/d' > {output.norm}
           
           
           """

rule somatic_merge_chrom:
    input: vcf=expand(os.path.join(output_somatic_snpindels,"{{vc_out}}","chrom_split","{{samples}}.{chroms}.vcf"), chroms=chroms),
    output: vcf=os.path.join(output_somatic_snpindels,"{vc_out}","vcf","{samples}.collected.vcf"),
    params: tumorsample="{samples}",genome = config['references']['GENOME'],genomedict = config['references']['GENOMEDICT'],ver_gatk=config['tools']['gatk4']['version'],rname="merge"
    shell: """
           module load GATK/{params.ver_gatk}
           input_str="-I $(echo "{input.vcf}" | sed -e 's/ / -I /g')"
           gatk --java-options "-Xmx30g" MergeVcfs -O "{output.vcf}" -D {params.genomedict} $input_str
           """


rule somatic_merge_callers:
    # input: vcf=expand(os.path.join(output_somatic_snpindels,"{vc_outdir}_out","vcf","{{samples}}.FINAL.vcf"), vc_outdir=caller_list)
    # output: mergedvcf=os.path.join(output_somatic_snpindels,"merged_somatic_variants","vcf","{samples}.FINAL.vcf"),
    input: vcf=expand(os.path.join(output_somatic_snpindels,"{vc_outdir}_out","vcf","{{samples}}.FINAL.norm.vcf"), vc_outdir=caller_list)
    output: mergedvcf=os.path.join(output_somatic_snpindels,"merged_somatic_variants","vcf","{samples}.FINAL.norm.vcf"),
    params: genome=config['references']['GENOME'],rodprioritylist=merge_callers_rodlist, variantsargs=lambda w: [merge_callers_args[w.samples]],ver_gatk=config['tools']['gatk3']['version'], rname="MergeSomaticCallers"
    shell: """
           if [ ! -d "$(dirname {output.mergedvcf})" ]; then
             mkdir -p "$(dirname {output.mergedvcf})"
           fi
           module load GATK/{params.ver_gatk}
           input_str="--variant $(echo "{input.vcf}" | sed -e 's/ / --variant /g')"
           nthread="$((SLURM_CPUS_PER_TASK-1))"
           java -Xmx60g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T CombineVariants -R {params.genome} -nt $nthread --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED --genotypemergeoption PRIORITIZE --rod_priority_list {params.rodprioritylist} --minimumN 1 -o {output.mergedvcf} {params.variantsargs}
           """
           
rule somatic_mafs:
  # input: filtered_vcf=os.path.join(output_somatic_snpindels,"{vc_outdir}","vcf","{samples}.FINAL.vcf"),
  input: filtered_vcf=os.path.join(output_somatic_snpindels,"{vc_outdir}","vcf","{samples}.FINAL.norm.vcf"),
         vcf2maf_script=VCF2MAF_WRAPPER
  output: maf=os.path.join(output_somatic_snpindels,"{vc_outdir}","maf","{samples}.maf")
  params: tumorsample="{samples}",genome=config['references']['MAF_GENOME'],filtervcf=config['references']['MAF_FILTERVCF'],rname="pl:vcf2maf"
  shell:
    """
    echo "Converting to MAF..."
    bash {input.vcf2maf_script} --vcf {input.filtered_vcf} --maf {output.maf} --tid {params.tumorsample} --genome {params.genome} --threads "$((SLURM_CPUS_PER_TASK-1))" --info "set"
    echo "Done converting to MAF..."
    date
    """
  
rule collect_cohort_mafs:
  input: mafs=expand(os.path.join(output_somatic_snpindels,"{{vc_outdir}}","maf","{samples}"+".maf"), samples=samples_for_caller_merge)
  output: maf=os.path.join(output_somatic_snpindels,"{vc_outdir}","cohort_summary","all_somatic_variants.maf")
  params: rname="combine_maf"
  shell:"""    
    date
    echo "Combining MAFs..."
    awk 'FNR==1 && NR!=1 {{ while (/^#/||/^Hugo_Symbol/) getline; }} 1 {{print}}' {input.mafs} > {output.maf}
    echo "Done..."
    date
  """
  