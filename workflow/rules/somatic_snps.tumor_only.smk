# Somatic SNP calling rules for tumor only samples
rule mutect2_single:
    input:
        tumor = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam")
    output:
        vcf = os.path.join(output_somatic_snpindels, "mutect2_out", "chrom_split", "{samples}.{chroms}.vcf"),
        read_orientation_file = os.path.join(output_somatic_snpindels, "mutect2_out", "chrom_split", "{samples}.{chroms}.f1r2.tar.gz"),
        statsfiles = os.path.join(output_somatic_snpindels, "mutect2_out", "chrom_split", "{samples}.{chroms}.vcf.stats")
    params:
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        pon = config['references']['PON'],
        germsource = config['references']['GNOMAD'],
        ver_gatk = config['tools']['gatk4']['version'],
        rname = 'mutect2'
    threads: 2
    envmodules:
        'GATK/4.2.0.0'
    container:
        config['images']['wes_base']
    shell: """
    if [ ! -d "$(dirname {output.vcf})" ]; then
        mkdir -p "$(dirname {output.vcf})"
    fi

    gatk Mutect2 \\
        -R {params.genome} \\
        -I {input.tumor} \\
        --panel-of-normals {params.pon} \\
        --germline-resource {params.germsource} \\
        -L {wildcards.chroms} \\
        -O {output.vcf} \\
        --f1r2-tar-gz {output.read_orientation_file} \\
        --independent-mates
    """


rule pileup_single:
    input:
        tumor = os.path.join(output_bamdir, "final_bams", "{samples}.bam"),
        bai = os.path.join(output_bamdir, "final_bams", "{samples}.bai"),
        intervals = intervals_file
    output:
        pileup = temp(os.path.join(output_somatic_snpindels, "mutect2_out", "pileup_summaries", "{samples}.pileup.table")),
    params:
        genome = config['references']['GENOME'],
        germsource = config['references']['1000GSNP'],
        ver_gatk = config['tools']['gatk4']['version'],
        chroms = chroms,
        rname = 'pileup',
        tmpdir = '/lscratch/$SLURM_JOBID'
    envmodules:
        'GATK/4.2.0.0'
    container:
        config['images']['wes_base']
    shell: """
    gatk --java-options "-Xmx10g -Djava.io.tmpdir={params.tmpdir}" GetPileupSummaries \\
        -R {params.genome} \\
        -I {input.tumor} \\
        -V {params.germsource} \\
        -L {input.intervals} \\
        -O {output.pileup}
    """


rule contamination_single:
    input: 
        pileup = os.path.join(output_somatic_snpindels, "mutect2_out", "pileup_summaries", "{samples}.pileup.table")
    output:
        tumor_summary = os.path.join(output_somatic_base, "qc", "gatk_contamination", "{samples}.contamination.table")
    params:
        genome = config['references']['GENOME'],
        germsource = config['references']['1000GSNP'],
        ver_gatk = config['tools']['gatk4']['version'],
        chroms = chroms, 
        rname = 'contamination'
    envmodules:
        'GATK/4.2.0.0'
    container:
        config['images']['wes_base']
    shell: """
    gatk CalculateContamination \\
        -I {input.pileup} \\
        -O {output.tumor_summary}
    """


rule mutect_single:
    input:
        tumor = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam"),
    output:
        vcf = os.path.join(output_somatic_snpindels, "mutect_out", "chrom_split", "{samples}.{chroms}.vcf"),
        stats = os.path.join(output_somatic_snpindels, "mutect_out", "chrom_split", "{samples}.{chroms}.stats.out"),
    params:
        genome = config['references']['GENOME'],
        pon = config['references']['PON'],
        cosmic = config['references']['COSMIC'],
        dbsnp = config['references']['DBSNP'],
        ver_mutect = config['tools']['mutect']['version'],
        rname = 'mutect',
        tmpdir = '/lscratch/$SLURM_JOBID'
    envmodules:
        'muTect/1.1.7'
    container:
        config['images']['mutect']
    shell: """
    if [ ! -d "$(dirname {output.vcf})" ]; then 
        mkdir -p "$(dirname {output.vcf})"
    fi

    java -Xmx8g -Djava.io.tmpdir={params.tmpdir} -jar ${{MUTECT_JAR}} \\
        --analysis_type MuTect \\
        --reference_sequence {params.genome} \\
        --normal_panel {params.pon} \\
        --vcf {output.vcf} \\
        --cosmic {params.cosmic} \\
        --dbsnp {params.dbsnp} \\
        -L {wildcards.chroms} \\
        --disable_auto_index_creation_and_locking_when_reading_rods \\
        --input_file:tumor {input.tumor} \\
        --out {output.stats} \\
        -rf BadCigar
    """


rule mutect_filter_single:
    input:
        vcf = os.path.join(output_somatic_snpindels, "mutect_out", "vcf", "{samples}.collected.vcf"),
    output:
        final = os.path.join(output_somatic_snpindels, "mutect_out", "vcf", "{samples}.FINAL.vcf"),
        norm = os.path.join(output_somatic_snpindels, "mutect_out", "vcf", "{samples}.FINAL.norm.vcf"),
    params:
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        ver_gatk = config['tools']['gatk4']['version'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'mutect_filter',
        tmpdir = '/lscratch/$SLURM_JOBID',
    threads: 4
    envmodules:
        'GATK/4.2.0.0',
        'bcftools/1.9'
    container:
        config['images']['wes_base']
    shell: """
    gatk SelectVariants \\
        -R {params.genome} \\
        --variant {input.vcf} \\
        --exclude-filtered \\
        --output {output.final}

    # VarScan can output ambiguous IUPAC bases/codes
    # the awk one-liner resets them to N, from:
    # https://github.com/fpbarthel/GLASS/issues/23
    bcftools sort -T {params.tmpdir} "{output.final}" \\
        | bcftools norm --threads {threads} --check-ref s -f {params.genome} -O v \\
        | awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' \\
        | sed '/^$/d' > {output.norm}
    """


rule vardict_single:
    input: 
        tumor = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam"),
    output:
        vcf = os.path.join(output_somatic_snpindels, "vardict_out", "chrom_split", "{samples}.{chroms}.vcf"),
    params: 
        genome = config['references']['GENOME'],
        targets = exome_targets_bed,
        pon = config['references']['PON'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'vardict'
    envmodules:
        'R/3.6.1',
        'samtools/1.8'
    container:
        config['images']['wes_base']
    shell: """
    if [ ! -d "$(dirname {output.vcf})" ]; then 
        mkdir -p "$(dirname {output.vcf})"
    fi

    VarDict \\
        -G {params.genome} \\
        -f 0.05 \\
        -x 500 \\
        --nosv \\
        -b {input.tumor} \\
        -t \\
        -Q 20 \\
        -c 1 \\
        -S 2 \\
        -E 3 {params.targets} \\
        | teststrandbias.R \\
        | var2vcf_valid.pl \\
            -N {wildcards.samples} \\
            -Q 20 \\
            -d 10 \\
            -v 6 \\
            -S \\
            -E \\
            -f 0.05 > {output.vcf}
    """


rule vardict_filter_single:
    input: 
        vcf = os.path.join(output_somatic_snpindels, "vardict_out", "vcf", "{samples}.collected.vcf"),
    output:
        final = os.path.join(output_somatic_snpindels, "vardict_out", "vcf", "{samples}.FINAL.vcf"),
        norm = os.path.join(output_somatic_snpindels, "vardict_out", "vcf", "{samples}.FINAL.norm.vcf"),
    params:
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        targets = exome_targets_bed,
        pon = config['references']['PON'],
        ver_gatk = config['tools']['gatk4']['version'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'vardict_filter', 
        tmpdir = '/lscratch/$SLURM_JOBID'
    threads: 4
    envmodules:
        'bcftools/1.9',
        'GATK/4.2.0.0'
    container:
        config['images']['wes_base']
    shell: """
    gatk SelectVariants \\
        -R {params.genome} \\
        --variant {input.vcf} \\
        --discordance {params.pon} \\
        --exclude-filtered \\
        --output {output.final}
    
    # VarScan can output ambiguous IUPAC bases/codes
    # the awk one-liner resets them to N, from:
    # https://github.com/fpbarthel/GLASS/issues/23
    bcftools sort -T {params.tmpdir} "{output.final}" \\
        | bcftools norm --threads {threads} --check-ref s -f {params.genome} -O v \\
        | awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' \\
        | sed '/^$/d' > {output.norm}
    """


rule varscan_single:
    input:
        tumor = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam"),
    output:
        vcf = os.path.join(output_somatic_snpindels, "varscan_out", "chrom_split", "{samples}.{chroms}.vcf"),
    params:
        genome = config['references']['GENOME'],
        ver_varscan = config['tools']['varscan']['version'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname='varscan'
    threads: 4
    envmodules:
        'VarScan/2.4.3',
        'bcftools/1.9'
    container:
        config['images']['wes_base']
    shell: """
    if [ ! -d "$(dirname {output.vcf})" ]; then
        mkdir -p "$(dirname {output.vcf})"
    fi

    varscan_opts="--strand-filter 0 --min-var-freq 0.01 --output-vcf 1 --variants 1"
    pileup_cmd="samtools mpileup -d 100000 -q 15 -Q 15 -f {params.genome} {input.tumor}"
    varscan_cmd="varscan mpileup2cns <($pileup_cmd) $varscan_opts"
    eval "$varscan_cmd > {output.vcf}.gz"
    
    # Create a vanilla empty VCF if no variants found
    if [ ! -s {output.vcf}.gz ]; then
        echo -e "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1" | bgzip -c > {output.vcf}.gz
    fi
    eval "bcftools view -U {output.vcf}.gz > {output.vcf}"
    """


rule varscan_filter_single:
    input:
        vcf = os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.collected.vcf"),
    output:
        filtered = temp(os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.filtered.vcf")),
        filtered1 = temp(os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.filtered.1.vcf")),
        samplesfile = temp(os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.FINAL.vcf.samples")),
        final = os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.FINAL.vcf"),
        norm = os.path.join(output_somatic_snpindels, "varscan_out", "vcf", "{samples}.FINAL.norm.vcf"),
    params:
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        pon = config['references']['PON'],
        basedir = BASEDIR,
        filter_settings = config['tools']['varscan']['filter_settings'],
        ver_gatk = config['tools']['gatk4']['version'],
        ver_varscan = config['tools']['varscan']['version'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'varscan_filter',
        tmpdir = '/lscratch/$SLURM_JOBID',
    threads: 4
    envmodules:
        'VarScan/2.4.3',
        'GATK/4.2.0.0',
        'bcftools/1.9'
    container:
        config['images']['wes_base']
    shell: """
    varscan filter \\
        {input.vcf} \\
        {params.filter_settings} > {output.filtered1}

    gatk SelectVariants \\
        -R {params.genome} \\
        --variant {output.filtered1} \\
        --discordance {params.pon} \\
        --exclude-filtered \\
        --output {output.filtered}

    samplesFile="{output.samplesfile}"
    echo -e "TUMOR\t{params.tumorsample}\n" > "{output.samplesfile}"

    bcftools reheader \\
        -o "{output.final}" \\
        -s "{output.samplesfile}" \\
        "{output.filtered}"

    # VarScan can output ambiguous IUPAC bases/codes
    # the awk one-liner resets them to N, from:
    # https://github.com/fpbarthel/GLASS/issues/23
    bcftools sort -T {params.tmpdir} "{output.final}" \\
        | bcftools norm --threads {threads} --check-ref s -f {params.genome} -O v \\
        | awk '{{gsub(/\y[W|K|Y|R|S|M]\y/,"N",$4); OFS = "\t"; print}}' \\
        | sed '/^$/d' > {output.norm}
    """