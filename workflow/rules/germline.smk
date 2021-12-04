# Rules for germline variant calling
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
        bam = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
        bai = os.path.join(output_bamdir,"final_bams","{samples}.bai"),
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
    message: "Running GATK4 HaplotypeCaller on '{input.bam}' input file"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base']
    shell:
        """
        myoutdir="$(dirname {output.gzvcf})"
        if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
         
        gatk --java-options '-Xmx24g' HaplotypeCaller \\
            --reference {params.genome} \\
            --input {input.bam} \\
            --use-jdk-inflater \\
            --use-jdk-deflater \\
            --emit-ref-confidence GVCF \\
            --annotation-group StandardAnnotation \\
            --annotation-group AS_StandardAnnotation \\
            --dbsnp {params.snpsites} \\
            --output {output.gzvcf} \\
            --intervals {params.chrom} \\
            --max-alternate-alleles 3
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
    params: 
        genome = config['references']['GENOME'],
        ver_gatk=config['tools']['gatk4']['version'],
        rname = "mergegvcfs"
    message: "Running GATK4 CombineGVCFs on '{input.gzvcf}' input file"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base'] 
    shell:
        """
        input_str="--variant $(echo "{input.gzvcf}" | sed -e 's/ / --variant /g')"
        
        gatk --java-options '-Xmx24g' CombineGVCFs \\
            --reference {params.genome} \\
            --annotation-group StandardAnnotation \\
            --annotation-group AS_StandardAnnotation \\
            $input_str \\
            --output {output.gzvcf} \\
            --intervals {wildcards.chroms} \\
            --use-jdk-inflater \\
            --use-jdk-deflater
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
    params:
        genome = config['references']['GENOME'],
        snpsites=config['references']['DBSNP'],
        chr="{chroms}",
        ver_gatk=config['tools']['gatk4']['version'],
        rname = "genotype"
    message: "Running GATK4 GenotypeGVCFs on '{input.gzvcf}' input file"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base']
    shell:
        """
        myoutdir="$(dirname {output.vcf})"
        if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
        
        gatk --java-options '-Xmx96g' GenotypeGVCFs \\
            --reference {params.genome} \\
            --use-jdk-inflater \\
            --use-jdk-deflater \\
            --annotation-group StandardAnnotation \\
            --annotation-group AS_StandardAnnotation \\
            --dbsnp {params.snpsites} \\
            --output {output.vcf} \\
            --variant {input.gzvcf} \\
            --intervals {params.chr}
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
        clist = os.path.join(output_germline_base,"VCF","by_chrom","raw_variants_byChrom.list"),
    params:
        rname = "merge_chrom", genome = config['references']['GENOME']
    message: "Running GATK4 MergeVcfs on all chrom split VCF files"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base'] 
    shell:
        """
        # Avoids ARG_MAX issue which limits max length of a command
        ls --color=never -d $(dirname "{output.clist}")/raw_variants.*.vcf.gz > "{output.clist}"

        gatk MergeVcfs \\
            -R {params.genome} \\
            --INPUT {output.clist} \\
            --OUTPUT {output.vcf}
        """


rule vqsr_snp:
    """
    Variant Quality Score Recalibration (VQSR) for SNPs
    @Input:
        Multi-sample gVCF
    @Output:
        Calculated info for applying SNP recalibration
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz")
    output:
        recal = os.path.join(output_germline_base,"vqsr","SNP.output.AS.recal"),
        tranches = os.path.join(output_germline_base,"vqsr","SNP.output.AS.tranches"),
        rscript = os.path.join(output_germline_base,"vqsr","SNP.output.plots.AS.R")
    params: 
        rname = "vqsr_snp",
        genome=config['references']['GENOME'],
        dbsnp=config['references']['DBSNP'],
        onekg=config['references']['1000GSNP'],
        hapmap=config['references']['HAPMAP'],
        omni=config['references']['OMNI'],
        gaussians = config['references']['VQSR']['SNP']['MAX_GAUSSIANS']
    message: "Running GATK4 VQSR for germline SNPs"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base'] 
    shell:
        """
        gatk --java-options '-Xmx24g' VariantRecalibrator \\
            -V {input.vcf} \\
            -O {output.recal} --tranches-file {output.tranches} --rscript-file {output.rscript} \\
            -mode SNP --trust-all-polymorphic \\
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \\
            --max-gaussians {params.gaussians} --reference {params.genome} --use-jdk-inflater --use-jdk-deflater -AS \\
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \\
            --resource:omni,known=false,training=true,truth=false,prior=12.0 {params.omni} \\
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.onekg} \\
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \\
            -an QD -an DP -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR 
        """

rule vqsr_indel:
    """
    Variant Quality Score Recalibration (VQSR) for INDELs
    @Input:
        Multi-sample gVCF
    @Output:
        Calculated info for applying INDEL recalibration
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
    output:
        recal = os.path.join(output_germline_base,"vqsr","INDEL.output.AS.recal"),
        tranches = os.path.join(output_germline_base,"vqsr","INDEL.output.AS.tranches"),
        rscript = os.path.join(output_germline_base,"vqsr","INDEL.output.plots.AS.R")
    params: 
        rname = "vqsr_indel",
        genome = config['references']['GENOME'],
        mills=config['references']['MILLS'],
        dbsnp=config['references']['DBSNP'],
        axiom=config['references']['AXIOM'],
        gaussians = config['references']['VQSR']['INDEL']['MAX_GAUSSIANS']
    message: "Running GATK4 VQSR for germline INDELs"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base'] 
    shell:
        """
        gatk --java-options '-Xmx24g' VariantRecalibrator \\
            -V {input.vcf} \\
            -mode INDEL --trust-all-polymorphic \\
            -O {output.recal} --tranches-file {output.tranches} --rscript-file {output.rscript} \\
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \\
            --reference {params.genome} --use-jdk-inflater --use-jdk-deflater -AS \\
            --resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mills} \\
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \\
            --resource:axiomPoly,known=false,training=true,truth=false,prior=10 {params.axiom} \\
            -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum --max-gaussians {params.gaussians}
        """

rule apply_snp:
    """
    Apply Variant Quality Score Recalibration (VQSR) for SNPs
    @Input:
        Multi-sample gVCF
        Recalibration file from GATK4 VariantRecalibrator in SNP mode
        Tranches file from GATK4 VariantRecalibrator in SNP mode
    @Output:
        gVCF with VQSR annotations
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
        recal = os.path.join(output_germline_base,"vqsr","SNP.output.AS.recal"),
        tranches = os.path.join(output_germline_base,"vqsr","SNP.output.AS.tranches",)
    output:
        vcf = os.path.join(output_germline_base,"vqsr","snps_recal_variants.vcf.gz"),
    params: 
        rname = "apply_snps",
        genome = config['references']['GENOME'],
        truth_level = config['references']['VQSR']['SNP']['TRUTH_SENSITIVITY']
    message: "Applying VQSR for germline SNPs"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base'] 
    shell:
        """
        gatk --java-options '-Xmx24g' ApplyVQSR \\
            -V {input.vcf} -mode SNP \\
            --recal-file {input.recal} --tranches-file {input.tranches} \\
            -O {output.vcf} \\
            --truth-sensitivity-filter-level {params.truth_level} \\
            --create-output-variant-index true \\
            --reference {params.genome} -AS --use-jdk-inflater --use-jdk-deflater
        """

rule apply_indel:
    """
    Apply Variant Quality Score Recalibration (VQSR) for INDELs
    @Input:
        Multi-sample gVCF
        Recalibration file from GATK4 VariantRecalibrator in INDEL mode
        Tranches file from GATK4 VariantRecalibrator in INDEL mode
    @Output:
        gVCF with VQSR annotations
    """
    input: 
        vcf = os.path.join(output_germline_base,"vqsr","snps_recal_variants.vcf.gz"),
        recal = os.path.join(output_germline_base,"vqsr","INDEL.output.AS.recal"),
        tranches = os.path.join(output_germline_base,"vqsr","INDEL.output.AS.tranches"),
    output:
        vcf = os.path.join(output_germline_base,"VCF","vqsr_variants.vcf.gz"),
    params: 
        rname = "apply_indels",
        genome = config['references']['GENOME'],
        truth_level = config['references']['VQSR']['INDEL']['TRUTH_SENSITIVITY']
    message: "Applying VQSR for germline INDELs"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base'] 
    shell:
        """
        gatk --java-options '-Xmx24g' ApplyVQSR \\
            -V {input.vcf} -mode INDEL \\
            --recal-file {input.recal} --tranches-file {input.tranches} \\
            -O {output.vcf} \\
            --truth-sensitivity-filter-level {params.truth_level} \\
            --reference {params.genome} -AS --use-jdk-inflater --use-jdk-deflater
        """

rule gtype_refinement:
    """
    Genotype refinement for germline variants
    @Input:
        Multi-sample gVCF
        Recalibration file from GATK4 VariantRecalibrator in INDEL mode
        Tranches file from GATK4 VariantRecalibrator in INDEL mode
    @Output:
        gVCF with VQSR annotations
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","vqsr_variants.vcf.gz"),
    output:
        vcf = temp(os.path.join(output_germline_base,"vqsr","snps_and_indels_recal_refinement_variants.vcf.gz")),
        gtfix_vcf = os.path.join(output_germline_base,"VCF","refined_germline_variants.vcf.gz"),
    params: 
        rname = "gtype_refinement",genome = config['references']['GENOME'],onekg = config['references']['1000GSNP']
    shell:
        """

        gatk --java-options '-Xmx24g' CalculateGenotypePosteriors \\
            -V {input.vcf} -supporting {params.onekg} \\
            -O {output.vcf} \\    
            --reference {params.genome} --use-jdk-inflater --use-jdk-deflater 

        bcftools +setGT {input.vcf} -O z -o {output.gtfix_vcf} -- -t a -n u
        tabix -p vcf {output.gtfix_vcf}
        """


rule Gatk_SelectVariants:
    """
    Make individual VCFs with variant sites only in that sample
    @Input:
        Multi-sample gVCF with all chromosomes combined
    @Output:
        Single-sample VCF with filtered germline variants
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","refined_germline_variants.vcf.gz"),
    output: 
        tmpvcf = temp(os.path.join(output_germline_base,"VCF","{samples}.select.vcf.gz")),
        vcf = os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz")
    params: 
        genome=config['references']['GENOME'], 
        Sname = "{samples}", 
        rname="varselect",
        ver_gatk=config['tools']['gatk4']['version'],
        targets=exome_targets_bed
    message: "Running GATK4 SelectVariants on '{input.vcf}' input file"
    envmodules: 'GATK/4.2.0.0'
    container: config['images']['wes_base']
    shell:
        """
        gatk SelectVariants \\
            -R {params.genome} \\
            --intervals {params.targets} \\
            --variant {input.vcf} \\
            --sample-name {params.Sname} \\
            --exclude-filtered \\
            --exclude-non-variants \\
            --keep-original-ac \\
            --output {output.tmpvcf}
        
        bcftools +setGT {output.tmpvcf} -O z -o {output.vcf} -- -t a -n u
        tabix -p vcf {output.vcf}
        """

