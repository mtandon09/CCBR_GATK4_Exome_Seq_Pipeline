
rule somalier_extract:
    input:  bam=os.path.join(output_bamdir,"final_bams","{samples}.bam"),
            bai=os.path.join(output_bamdir,"final_bams","{samples}.bai"),
    output: somalierOut=os.path.join(output_germline_base,"somalier","{samples}.somalier")
    params: ancestry_db=config['references']['SOMALIER']['ANCESTRY_DB'],somalier_container=config['references']['SOMALIER']['CONTAINER'],
            sites_vcf=config['references']['SOMALIER']['SITES_VCF'],genomeFasta=config['references']['GENOME'],rname="somalier_extract"
    # threads: 32
    shell:  """ 
        module load singularity  
	    ANCESTRY_DB={params.ancestry_db}
        CONTAINER={params.somalier_container}
        GENOME_FASTA={params.genomeFasta}
        SITES_VCF={params.sites_vcf}
	    bam={input.bam}
        OUT_DIR=$(dirname {output.somalierOut}) #Need to set this path somewhere
        #mkdir -p $OUT_DIR
        BAM_BASE=$(basename $bam)
        SOMALIER_CMDS="date; echo 'Extracting sites';
                    somalier extract -d /out/ --sites /refs/sites.vcf.gz -f /refs/genome.fa /mnt/$BAM_BASE;
                    date;"
            ## Run the commands in the image with binded paths to input/output files
            singularity exec --bind $GENOME_FASTA:/refs/genome.fa \
                                     --bind $GENOME_FASTA.fai:/refs/genome.fa.fai \
                                     --bind $SITES_VCF:/refs/sites.vcf.gz \
                                     --bind $SITES_VCF.tbi:/refs/sites.vcf.gz.tbi \
                                     --bind $ANCESTRY_DB:/refs/ancestry/ \
                                     --bind $OUT_DIR:/out \
                                     --bind $bam:/mnt/$BAM_BASE \
                                     --bind $bam.bai:/mnt/$BAM_BASE.bai \
                                     $CONTAINER \
                                     /bin/sh -c "$SOMALIER_CMDS"
            """
  
rule somalier_relatedness:
    input:
        somalier=expand(os.path.join(output_germline_base,"somalier","{samples}.somalier"), samples=samples),
    output:
        relatedness=os.path.join(output_germline_base,"somalier","relatedness.pairs.tsv"),
        relatednessSamples=os.path.join(output_germline_base,"somalier","relatedness.samples.tsv"),
        ancestry=os.path.join(output_germline_base,"somalier","ancestry.somalier-ancestry.tsv"),
    params: ancestry_db=config['references']['SOMALIER']['ANCESTRY_DB'],container=config['references']['SOMALIER']['CONTAINER'],
            sites_vcf=config['references']['SOMALIER']['SITES_VCF'],genomeFasta=config['references']['GENOME'],rname="somalier_relatedness"
    shell:  """
        module load singularity    
	    ANCESTRY_DB={params.ancestry_db}
        CONTAINER={params.container}
        GENOME_FASTA={params.genomeFasta}
        SITES_VCF={params.sites_vcf}
        OUT_DIR=$(dirname {output.relatedness}) #Need to set this path somewhere
        SOMALIER_CMDS="date;
                    somalier relate -o /out/relatedness /out/*.somalier;
                    somalier ancestry -o /out/ancestry --labels /refs/ancestry/ancestry-labels-1kg.tsv /refs/ancestry/*.somalier ++ /out/*.somalier;
                    echo 'Done';
                    date;"
            ## Run the commands in the image with binded paths to input/output files
            singularity exec --bind $GENOME_FASTA:/refs/genome.fa \
                                     --bind $GENOME_FASTA.fai:/refs/genome.fa.fai \
                                     --bind $SITES_VCF:/refs/sites.vcf.gz \
                                     --bind $SITES_VCF.tbi:/refs/sites.vcf.gz.tbi \
                                     --bind $ANCESTRY_DB:/refs/ancestry/ \
                                     --bind $OUT_DIR:/out \
                                     $CONTAINER \
                                     /bin/sh -c "$SOMALIER_CMDS"
            """

rule somalier_analysis:
    input:  somalierPairs=os.path.join(output_germline_base,"somalier","relatedness.pairs.tsv"),
            somalierSamples=os.path.join(output_germline_base,"somalier","relatedness.samples.tsv"),
            somalierAncestry=os.path.join(output_germline_base,"somalier","relatedness.pairs.tsv"),
    output: finalFileGender=os.path.join(output_germline_base,"predicted.genders.tsv"),
            finalFilePairs=os.path.join(output_germline_base,"predicted.pairs.tsv"),
            ancestoryPlot=os.path.join(output_germline_base,"sampleAncestryPCAPlot.html"),
            pairAncestoryHist=os.path.join(output_germline_base,"predictedPairsAncestry.pdf"),
    params: ver_R=config['tools']['R']['version'],
            script_path_gender=config['scripts']['genderPrediction'],
            script_path_samples=config['scripts']['combineSamples'],
            script_path_pca=config['scripts']['ancestry'],
            rname = "genderPrediction"
    shell: """
            module load R/{params.ver_R}
            
            Rscript {params.script_path_gender} {input.somalierSamples} {output.finalFileGender}
            
            Rscript {params.script_path_samples} {input.somalierPairs} {output.finalFilePairs}
            
            Rscript {params.script_path_pca} {input.somalierAncestry} {output.finalFilePairs} {output.ancestoryPlot} {output.pairAncestoryHist}
    """


# Quality-control related rules
rule fc_lane:
    """
    Quality-control step to get flowcell and lane information from FastQ file.
    FastQ files generated with older versions of Casava or downloaded from
    SRA have a different format than newer FastQ files generated with the
    current version of Casava. It is worth noting that FastQ files downloaded from SRA
    or FastQ files generated with Casava version < 1.8 do not have Flowcell
    IDs in its sequence indentifer. If a FastQ file does not have Flowcell IDs,
    the Machine or Instrument ID is grabbed instead.
    @Input:
        Raw FastQ R1 file (scatter)
    @Output:
        Text file containing information about the FastQ file
    """
    input:
        r1 = os.path.join(input_fqdir, "{samples}.R1.fastq.gz"),
    output:
        txt = os.path.join(output_fqdir,"{samples}.fastq.info.txt")
    params:
        rname = 'pl:fc_lane',
        get_flowcell_lanes = os.path.join("scripts", "get_flowcell_lanes.py"),
    shell: """
    module load python/2.7
    python {params.get_flowcell_lanes} \\
        {input.r1} \\
        {wildcards.samples} > {output.txt}
    """


rule fastq_screen:
    """
    Quality-control step to screen for different sources of contamination.
    FastQ Screen compares your sequencing data to a set of different reference
    genomes to determine if there is contamination. It allows a user to see if
    the composition of your library matches what you expect.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        FastQ Screen report and logfiles
    """
    input:
        fq1 = os.path.join(output_fqdir,"{samples}.R1.trimmed.fastq.gz"),
        fq2 = os.path.join(output_fqdir,"{samples}.R2.trimmed.fastq.gz")
    output:
        txt1 = os.path.join(output_qcdir,"FQscreen","{samples}.R1.trimmed_screen.txt"),
        txt2 = os.path.join(output_qcdir,"FQscreen","{samples}.R2.trimmed_screen.txt"),
        png1 = os.path.join(output_qcdir,"FQscreen","{samples}.R1.trimmed_screen.png"),
        png2 = os.path.join(output_qcdir,"FQscreen","{samples}.R2.trimmed_screen.png")
    params:
        rname  = "pl:fqscreen",
        outdir = os.path.join(output_qcdir,"FQscreen"),
        # Exposed Parameters: modify resources/fastq_screen.conf to change 
        # default locations to bowtie2 indices
        fastq_screen_config = config['references']['FASTQ_SCREEN_CONFIG'],
    threads: 24
    shell: """
    module load fastq_screen/0.14.1
    fastq_screen --conf {params.fastq_screen_config} \\
        --outdir {params.outdir} \\
        --threads {threads} \\
        --subset 1000000 \\
        --aligner bowtie2 \\
        --force \\
        {input.fq1} {input.fq2}
    """


rule kraken:
    """
    Quality-control step to assess for potential sources of microbial contamination.
    If there are high levels of microbial contamination, Kraken will provide an
    estimation of the taxonomic composition. Kraken is used in conjunction with
    Krona to produce an interactive reports.
    @Input:
        Trimmed FastQ files (scatter)
    @Output:
        Kraken logfile and interative krona report
    """
    input:
        fq1 = os.path.join(output_fqdir,"{samples}.R1.trimmed.fastq.gz"),
        fq2 = os.path.join(output_fqdir,"{samples}.R2.trimmed.fastq.gz")
    output:
        out  = os.path.join(output_qcdir,"kraken","{samples}.trimmed.kraken_bacteria.out.txt"),
        taxa = os.path.join(output_qcdir,"kraken","{samples}.trimmed.kraken_bacteria.taxa.txt"),
        html = os.path.join(output_qcdir,"kraken","{samples}.trimmed.kraken_bacteria.krona.html"),
    params:
        rname  ='pl:kraken',
        outdir = os.path.join(output_qcdir, "kraken"),
        bacdb  = config['references']['KRAKENBACDB'],
    threads: 24
    shell: """
    module load kraken/2.1.2
    module load kronatools/2.8
    # Copy kraken2 db to local node storage to reduce filesystem strain
    cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/
    kdb_base=$(basename {params.bacdb})
    kraken2 --db /lscratch/$SLURM_JOBID/${{kdb_base}} \\
        --threads {threads} --report {output.taxa} \\
        --output {output.out} \\
        --gzip-compressed \\
        --paired {input.fq1} {input.fq2}
    # Generate Krona Report
    cut -f2,3 {output.out} | \\
        ktImportTaxonomy - -o {output.html}
    """


rule fastqc_bam:
    """
    Quality-control step to assess sequencing quality of each sample.
    FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        FastQC report and zip file containing sequencing quality information
    """
    input:
        bam = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
    output:
        zipfile =  os.path.join(output_qcdir,"{samples}.recal_fastqc.zip"),
        report  =  os.path.join(output_qcdir,"{samples}.recal_fastqc.html")
    params:
        outdir = output_qcdir,
        rname  = "fastqc_bam",
    message: "Running FastQC with {threads} threads on '{input}' input file"
    threads: 8
    shell: """
    module load fastqc/0.11.9
    fastqc -t {threads} \\
        -f bam \\
        -o {params.outdir} \\
        {input.bam} 
    """


rule qualimap_bamqc:
    """
    Quality-control step to assess various post-alignment metrics 
    and a secondary method to calculate insert size. Please see
    QualiMap's website for more information about BAM QC:
    http://qualimap.conesalab.org/
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        Report containing post-aligment quality-control metrics
    """
    input:
        bam  = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
    output: 
        txt  = os.path.join(output_qcdir,"{samples}","genome_results.txt"),
        html = os.path.join(output_qcdir,"{samples}","qualimapReport.html")
    params:
        outdir = os.path.join(output_qcdir, "{samples}"),
        rname  = "qualibam"
    message: "Running QualiMap BAM QC with {threads} threads on '{input}' input file"
    threads: 8
    shell: """
    module load qualimap/2.2.1
    unset DISPLAY
    qualimap bamqc -bam {input.bam} \\
        --java-mem-size=48G \\
        -c -gd hg19 -ip \\
        -outdir {params.outdir} \\
        -outformat HTML \\
        -nt {threads} \\
        --skip-duplicated \\
        -nw 500 \\
        -p NON-STRAND-SPECIFIC
    """


rule samtools_flagstats:
    """
    Quality-control step to assess alignment quality. Flagstat provides 
    counts for each of 13 categories based primarily on bit flags in the 
    FLAG field. Information on the meaning of the flags is given in the 
    SAM specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    @Input:
        Recalibrated BAM file (scatter)
    @Output:
        Text file containing alignment statistics
    """
    input:
        bam  = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
    output:
        txt  = os.path.join(output_qcdir,"{samples}.samtools_flagstat.txt")
    params: 
        rname = "samtools_flagstats"
    message: "Running SAMtools flagstat on '{input}' input file"
    shell: """
    module load samtools/1.12
    samtools flagstat {input.bam} > {output.txt}
    """


rule vcftools:
    """
    Quality-control step to calculates a measure of heterozygosity on 
    a per-individual basis. The inbreeding coefficient, F, is estimated
    for each individual using a method of moments. Please see VCFtools
    documentation for more information: 
    https://vcftools.github.io/man_latest.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a measure of heterozygosity
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
    output: 
        het = os.path.join(output_qcdir,"raw_variants.het"),
    params: 
        prefix = os.path.join(output_qcdir,"raw_variants"),
        rname  = "vcftools",
    message: "Running VCFtools on '{input.vcf}' input file"
    shell: """
    module load vcftools/0.1.16
    vcftools --gzvcf {input.vcf} --het --out {params.prefix}
    """


rule collectvariantcallmetrics:
    """
    Quality-control step to collect summary metrics about snps and indels
    called in a multisample VCF file. Please see the Broad's documentation
    for more information about each field in the generated log file:
    https://broadinstitute.github.io/picard/picard-metric-definitions.html
    @Input:
        Multi-sample gVCF file (indirect-gather-due-to-aggregation)
    @Output:
        Text file containing a collection of metrics relating to snps and indels 
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","raw_variants.vcf.gz"),
    output: 
        metrics = os.path.join(output_qcdir,"raw_variants.variant_calling_detail_metrics"),
    params: 
        dbsnp=config['references']['DBSNP'],
        prefix = os.path.join(output_qcdir,"raw_variants"),
        rname="varcallmetrics",
    message: "Running Picard CollectVariantCallingMetrics on '{input.vcf}' input file"
    shell: """
    module load picard/2.20.8
    java -Xmx24g -jar ${{PICARDJARPATH}}/picard.jar \\
        CollectVariantCallingMetrics \\
        INPUT={input.vcf} \\
        OUTPUT={params.prefix} \\
        DBSNP={params.dbsnp} Validation_Stringency=SILENT
    """


rule bcftools_stats:
    """
    Quality-control step to collect summary statistics from bcftools stats.
    When bcftools stats is run with one VCF file then stats by non-reference
    allele frequency, depth distribution, stats by quality and per-sample 
    counts, singleton statsistics are calculated. Please see bcftools' 
    documentation for more information: 
    http://samtools.github.io/bcftools/bcftools.html#stats
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Text file containing a collection of summary statistics
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz"),
    output: 
        txt = os.path.join(output_qcdir,"{samples}.germline.bcftools_stats.txt"),
    params: 
        rname="bcfstats",
    shell: """
    module load bcftools/1.9
    bcftools stats {input.vcf} > {output.txt}
    """

rule gatk_varianteval:
    """
    Quality-control step to calculate various quality control metrics from a 
    variant callset. These metrics include the number of raw or filtered SNP 
    counts; ratio of transition mutations to transversions; concordance of a
    particular sample's calls to a genotyping chip; number of s per sample.
    Please see GATK's documentation for more information: 
    https://gatk.broadinstitute.org/hc/en-us/articles/360040507171-VariantEval
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Evaluation table containing a collection of summary statistics
    """
    input: 
        vcf = os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz"), 
    output: 
        grp = os.path.join(output_qcdir,"{samples}.germline.eval.grp"),
    params:
        rname    = "vareval",
        genome   = config['references']['GENOME'],
        dbsnp    = config['references']['DBSNP'],
        ver_gatk = config['tools']['gatk4']['version']
    threads: 16
    shell: """
    module load GATK/{params.ver_gatk}
    gatk --java-options '-Xmx12g -XX:ParallelGCThreads={threads}' VariantEval \\
        -R {params.genome} \\
        -O {output.grp} \\
        --dbsnp {params.dbsnp} \\
        --eval {input.vcf} 
    """


rule snpeff:
    """
    Data processing and quality-control step to annotate variants, predict its
    functional effects, and collect various summary statistics about variants and
    their annotations. Please see SnpEff's documentation for more information: 
    https://pcingola.github.io/SnpEff/
    @Input:
        Per sample gVCF file (scatter)
    @Output:
        Evaluation table containing a collection of summary statistics
    """
    input:  
        vcf = os.path.join(output_germline_base,"VCF","{samples}.germline.vcf.gz")
    output: 
        vcf  = os.path.join(output_qcdir,"{samples}.germline.snpeff.ann.vcf"),
        csv  = os.path.join(output_qcdir,"{samples}.germline.snpeff.ann.csv"),
        html = os.path.join(output_qcdir,"{samples}.germline.snpeff.ann.html"),
    params: 
        rname  = "snpeff",
        genome = config['references']['SNPEFF_GENOME'],
        config = config['references']['SNPEFF_CONFIG']
    shell: """
    module load snpEff/4.3t
    java -Xmx12g -jar $SNPEFF_JAR \\
        -v -canon -c {params.config} \\
        -csvstats {output.csv} \\
        -stats {output.html} \\
        {params.genome} \\
        {input.vcf} > {output.vcf}
    """


rule multiqc:
    """
    Reporting step to aggregate sample summary statistics and quality-control
    information across all samples. This will be one of the last steps of the 
    pipeline. The inputs listed here are to ensure that this step runs last. 
    During runtime, MultiQC will recurively crawl through the working directory
    and parse files that it supports.
    @Input:
        List of files to ensure this step runs last (gather)
    @Output:
        Interactive MulitQC report and a QC metadata table
    """
    input:  
        expand(os.path.join(output_fqdir,"{samples}.fastq.info.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"FQscreen","{samples}.R2.trimmed_screen.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"kraken","{samples}.trimmed.kraken_bacteria.krona.html"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.recal_fastqc.zip"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}","genome_results.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.samtools_flagstat.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.germline.bcftools_stats.txt"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.germline.eval.grp"), samples=samples),
        expand(os.path.join(output_qcdir,"{samples}.germline.snpeff.ann.html"), samples=samples),
        os.path.join(output_qcdir,"raw_variants.het"), 
        os.path.join(output_qcdir,"raw_variants.variant_calling_detail_metrics"),
        os.path.join(output_germline_base,"somalier","ancestry.somalier-ancestry.tsv"),
    output: 
        report  = os.path.join(output_qcdir,"MultiQC_Report.html"),
    params: 
        rname  = "multiqc",
        workdir = os.path.join(BASEDIR)
    shell: """
    module load multiqc/1.11
    multiqc --ignore '*/.singularity/*' \\
        --ignore '*/*/*/*/*/*/*/*/pyflow.data/*' \\
        -f --interactive \\
        -n {output.report} \\
        {params.workdir}
    """