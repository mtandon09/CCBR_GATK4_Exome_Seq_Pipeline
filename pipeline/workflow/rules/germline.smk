
rule haplotypecaller:
    """
    Germline variant calling. This can be done independently across the
    genome, so we're splitting it up by chromosome.
    @Input:
        Aligned reads in BAM format
    @Output:
        Single-sample gVCF
    """
    input: 
        bam=os.path.join(output_bamdir,"final_bams","{samples}.bam"),
        bai=os.path.join(output_bamdir,"final_bams","{samples}.bai"),
    output:
        gzvcf = temp(os.path.join(output_germline_base,"gVCFs","{samples}.{chroms}.g.vcf.gz")),
        index = temp(os.path.join(output_germline_base,"gVCFs","{samples}.{chroms}.g.vcf.gz.tbi")),
    params: 
        sample = "{samples}",
        genome = config['references']['GENOME'],
        snpsites=config['references']['DBSNP'],
        chrom="{chroms}",
        ver_gatk=config['tools']['gatk4']['version'],
        rname = "hapcaller"
    shell:
        """
        myoutdir="$(dirname {output.gzvcf})"
        if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
         
        module load GATK/{params.ver_gatk}
        gatk --java-options '-Xmx24g' HaplotypeCaller --reference {params.genome} --input {input.bam} --use-jdk-inflater --use-jdk-deflater --emit-ref-confidence GVCF --annotation-group StandardAnnotation --annotation-group AS_StandardAnnotation --dbsnp {params.snpsites} --output {output.gzvcf} --intervals {params.chrom} --max-alternate-alleles 3
        """

rule mergegvcfs:
    """
    Merge germline variants across samples
    @Input:
        Single-sample gVCFs, scattered across chromosomes
    @Output:
        Multi-sample gVCF, scattered across chromosomes
    """
    input: gzvcf = expand(os.path.join(output_germline_base,"gVCFs","{samples}.{{chroms}}.g.vcf.gz"),samples=samples),
           index = expand(os.path.join(output_germline_base,"gVCFs","{samples}.{{chroms}}.g.vcf.gz.tbi"),samples=samples),
           # list = "gVCFs/gVCFs.{chroms}.list",
    output:
        gzvcf = os.path.join(output_germline_base,"gVCFs","merged.{chroms}.g.vcf.gz"),
        index = os.path.join(output_germline_base,"gVCFs","merged.{chroms}.g.vcf.gz.tbi"),
    # group: "merge_and_genotype"
    params: 
        genome = config['references']['GENOME'],
        ver_gatk=config['tools']['gatk4']['version'],
        rname = "mergegvcfs"
    shell:
        """
        input_str="--variant $(echo "{input.gzvcf}" | sed -e 's/ / --variant /g')"
        
        module load GATK/{params.ver_gatk}
        gatk --java-options '-Xmx24g' CombineGVCFs --reference {params.genome} --annotation-group StandardAnnotation --annotation-group AS_StandardAnnotation $input_str --output {output.gzvcf} --intervals {wildcards.chroms} --use-jdk-inflater --use-jdk-deflater
        """

rule genotype:
    """
    Joint genotyping of germline variants
    @Input:
        Multi-sample gVCF, scattered across chromosomes
    @Output:
        Multi-sample gVCF, scattered across chromosomes (with joint genotyping updates)
    """
    input: 
        gzvcf = os.path.join(output_germline_base,"gVCFs","merged.{chroms}.g.vcf.gz"),
        index = os.path.join(output_germline_base,"gVCFs","merged.{chroms}.g.vcf.gz.tbi"),
    output:
        vcf = os.path.join(output_germline_base,"VCF","by_chrom","raw_variants.{chroms}.vcf.gz"),
    # group: "merge_and_genotype"
    params:
        genome = config['references']['GENOME'],
        snpsites=config['references']['DBSNP'],
        chr="{chroms}",
        ver_gatk=config['tools']['gatk4']['version'],
        rname = "genotype"
    shell:
        """
        myoutdir="$(dirname {output.vcf})"
        if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
        
        module load GATK/{params.ver_gatk}
        gatk --java-options '-Xmx96g' GenotypeGVCFs --reference {params.genome} --use-jdk-inflater --use-jdk-deflater --annotation-group StandardAnnotation --annotation-group AS_StandardAnnotation --dbsnp {params.snpsites} --output {output.vcf} --variant {input.gzvcf} --intervals {params.chr}
        """

rule germline_merge_chrom:
    """
    Combine joint genotyping from all chromosomes
    @Input:
        Multi-sample gVCF for all chromosomes
    @Output:
        Multi-sample gVCF with all chromosomes combined
    """
    input:
        expand(os.path.join(output_germline_base,"VCF","by_chrom","raw_variants.{chroms}.vcf.gz"), chroms=chroms),
    output:
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
        list = os.path.join(output_germline_base,"VCF","by_chrom","raw_variants_byChrom.list"),
    params:
        rname = "merge_chrom", genome = config['references']['GENOME']
    shell:
        """
        ls -d $(dirname "{output.list}")/raw_variants.*.vcf.gz > "{output.list}"
        module load GATK/4.1.7.0
        gatk MergeVcfs -R {params.genome} --INPUT="{output.list}" --OUTPUT={output.vcf}
        """

rule Gatk_SelectVariants:
    """
    Make individual VCFs with variant sites in that sample
    @Input:
        Multi-sample gVCF with all chromosomes combined
    @Output:
        Single-sample VCF with unfiltered germline variants
    """
	input: vcf=os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
	output: vcf=os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz")
	params: genome=config['references']['GENOME'], Sname = "{samples}", rname="varselect",ver_gatk=config['tools']['gatk4']['version'],targets=exome_targets_bed
	shell:"""
          module load GATK/{params.ver_gatk}
          gatk SelectVariants -R {params.genome} --intervals {params.targets} --variant {input.vcf} --sample-name {params.Sname} --exclude-filtered --exclude-non-variants --output {output.vcf}
	      """
          