# Common somatic SNP calling rules
rule split_bam_by_chrom:
    input:
        bam = os.path.join(output_bamdir, "final_bams", "{samples}.bam"),
        bai = os.path.join(output_bamdir, "final_bams", "{samples}.bam.bai"),
    output:
        split_bam = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bam"),
        split_bam_idx = os.path.join(output_bamdir, "chrom_split", "{samples}.{chroms}.split.bai")
    params:
        ver_samtools = config['tools']['samtools']['version'],
        rname='bam_split'
    threads: 4
    envmodules:
        'samtools/1.12'
    container:
       config['images']['wes_base'] 
    shell: """
    if [ ! -d "$(dirname {output.split_bam})" ]; then
      mkdir -p "$(dirname {output.split_bam})"
    fi

    samtools view \\
        -b \\
        -o {output.split_bam} \\
        -@ {threads} \\
        {input.bam} {wildcards.chroms}
    
    samtools index \\
        -@ {threads} \\
        {output.split_bam} {output.split_bam_idx}
    
    cp {output.split_bam_idx} {output.split_bam}.bai
    """


rule LearnReadOrientationModel:
    input:
        vcf = expand(os.path.join(output_somatic_snpindels, "mutect2_out", "chrom_split", "{{samples}}.{chroms}.vcf"), chroms=chroms),
        read_orientation_file = expand(os.path.join(output_somatic_snpindels, "mutect2_out", "chrom_split", "{{samples}}.{chroms}.f1r2.tar.gz"), chroms=chroms)
    output:
        model = os.path.join(output_somatic_snpindels, "mutect2_out", "read_orientation_data", "{samples}.read-orientation-model.tar.gz")
    params:
        genome = config['references']['GENOME'],
        ver_gatk = config['tools']['gatk4']['version'],
        rname = 'LearnReadOrientationModel'
    envmodules:
        'GATK/4.2.0.0'
    container:
        config['images']['wes_base']
    shell: """
    input_str="--input $(echo "{input.read_orientation_file}" | sed -e 's/ / --input /g')"
     
    gatk LearnReadOrientationModel \\
        --output {output.model} \\
        $input_str
    """


rule mutect2_filter:
    input:
        vcf = os.path.join(output_somatic_snpindels, "mutect2_out", "vcf", "{samples}.collected.vcf"),
        summary = os.path.join(output_somatic_base, "qc", "gatk_contamination", "{samples}.contamination.table"),
        model = os.path.join(output_somatic_snpindels, "mutect2_out", "read_orientation_data", "{samples}.read-orientation-model.tar.gz"),
        statsfiles = expand(os.path.join(output_somatic_snpindels, "mutect2_out", "chrom_split", "{{samples}}.{chroms}.vcf.stats"), chroms=chroms)
    output:
        marked_vcf = os.path.join(output_somatic_snpindels, "mutect2_out", "vcf", "{samples}.filtered.vcf"),
        final = os.path.join(output_somatic_snpindels, "mutect2_out", "vcf", "{samples}.FINAL.vcf"),
        norm = os.path.join(output_somatic_snpindels, "mutect2_out", "vcf", "{samples}.FINAL.norm.vcf"),
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        ver_gatk = config['tools']['gatk4']['version'],
        ver_bcftools = config['tools']['bcftools']['version'],
        rname = 'mutect2_filter',
        tmpdir = '/lscratch/$SLURM_JOBID',
    threads: 2
    envmodules:
        'GATK/4.2.0.0',
        'bcftools/1.9'
    container:
        config['images']['wes_base']
    shell: """
    statfiles="--stats $(echo "{input.statsfiles}" | sed -e 's/ / --stats /g')"
    
    gatk MergeMutectStats \\
        $statfiles \\
        -O {output.final}.stats

    gatk FilterMutectCalls \\
        -R {params.genome} \\
        -V {input.vcf} \\
        --ob-priors {input.model} \\
        --contamination-table {input.summary} \\
        -O {output.marked_vcf} \\
        --stats {output.final}.stats

    gatk SelectVariants \\
        -R {params.genome} \\
        --variant {output.marked_vcf} \\
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
           
    
rule somatic_merge_chrom:
    input:
        vcf = expand(os.path.join(output_somatic_snpindels, "{{vc_out}}", "chrom_split", "{{samples}}.{chroms}.vcf"), chroms=chroms),
    output:
        vcf = os.path.join(output_somatic_snpindels, "{vc_out}", "vcf", "{samples}.collected.vcf"),
    params:
        tumorsample = '{samples}',
        genome = config['references']['GENOME'],
        genomedict = config['references']['GENOMEDICT'],
        ver_gatk = config['tools']['gatk4']['version'],
        rname = 'merge'
    envmodules:
        'GATK/4.2.0.0'
    container:
        config['images']['wes_base']
    shell: """
    input_str="-I $(echo "{input.vcf}" | sed -e 's/ / -I /g')"

    gatk --java-options "-Xmx30g" MergeVcfs \\
        -O "{output.vcf}" \\
        -D {params.genomedict} \\
        $input_str
    """


rule somatic_merge_callers:
    input:
        vcf = expand(os.path.join(output_somatic_snpindels, "{vc_outdir}_out", "vcf", "{{samples}}.FINAL.norm.vcf"), vc_outdir=caller_list)
    output: 
        mergedvcf = os.path.join(output_somatic_snpindels, "merged_somatic_variants", "vcf", "{samples}.FINAL.norm.vcf"),
    params:
        genome = config['references']['GENOME'],
        rodprioritylist = merge_callers_rodlist,
        variantsargs = lambda w: [merge_callers_args[w.samples]],
        ver_gatk = config['tools']['gatk3']['version'],
        rname = 'MergeSomaticCallers',
        tmpdir = '/lscratch/$SLURM_JOBID'
    threads: 4
    envmodules:
        'GATK/3.8-1'
    container:
        config['images']['wes_base']
    shell: """
    if [ ! -d "$(dirname {output.mergedvcf})" ]; then
      mkdir -p "$(dirname {output.mergedvcf})"
    fi

    input_str="--variant $(echo "{input.vcf}" | sed -e 's/ / --variant /g')"

    java -Xmx60g -Djava.io.tmpdir={params.tmpdir} -jar $GATK_JAR -T CombineVariants \\
        -R {params.genome} \\
        -nt {threads} \\
        --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \\
        --genotypemergeoption PRIORITIZE \\
        --rod_priority_list {params.rodprioritylist} \\
        --minimumN 1 \\
        -o {output.mergedvcf} \\
        {params.variantsargs}
    """


rule somatic_mafs:
    input:
        filtered_vcf = os.path.join(output_somatic_snpindels, "{vc_outdir}", "vcf", "{samples}.FINAL.norm.vcf")
    output:
        maf = os.path.join(output_somatic_snpindels, "{vc_outdir}", "maf", "{samples}.maf")
    params: 
        tumorsample = '{samples}',
        genome = config['references']['MAF_GENOME'],
        filtervcf = config['references']['MAF_FILTERVCF'],
        bundle = config['references']['VCF2MAF']['VEPRESOURCEBUNDLEPATH'],
        rname = 'vcf2maf',
        vcf2maf_script = VCF2MAF_WRAPPER
    threads: 4
    # container:
    #     config['images']['vcf2maf'] 
    shell: """
    echo "Converting to MAF..."
      
    bash {params.vcf2maf_script} \\
        --vcf {input.filtered_vcf} \\
        --maf {output.maf} \\
        --tid {params.tumorsample} \\
        --genome {params.genome} \\
        --threads {threads} \\
        --vepresourcebundlepath {params.bundle} \\
        --info "set"
    echo "Done converting to MAF..."
    """


rule collect_cohort_mafs:
    input: 
        mafs = expand(os.path.join(output_somatic_snpindels, "{{vc_outdir}}", "maf", "{samples}"+".maf"), samples=samples_for_caller_merge)
    output: 
        maf = os.path.join(output_somatic_snpindels, "{vc_outdir}", "cohort_summary", "all_somatic_variants.maf")
    params:
        rname = 'combine_maf'
    shell: """
    echo "Combining MAFs..."
    head -2 {input.mafs[0]} > {output.maf}
    awk 'FNR>2 {{print}}' {input.mafs} >> {output.maf}
    """