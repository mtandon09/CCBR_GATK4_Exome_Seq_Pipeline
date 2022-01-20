# Rules for structural variant calling
rule manta_somatic:
    input: 
        normal = lambda w: [os.path.join(output_bamdir,"final_bams", pairs_dict[w.samples] + ".bam")],
        tumor = os.path.join(output_bamdir,"final_bams","{samples}.bam"),
    output: 
        sv_vcf=os.path.join(output_somatic_sv, "manta_out","{samples}","results/variants/somaticSV.vcf.gz")
    params:
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = "{samples}",
        genome = config['references']['GENOME'],
        outdir = output_somatic_sv,
        rname = 'manta_somatic',
    envmodules: 
        'manta/1.6.0',
    container: config['images']['wes_base']
    shell: """
    myoutdir="{params.outdir}/manta_out/{params.tumorsample}"
    if [ -d "$myoutdir" ]; then rm -r "$myoutdir"; fi
    mkdir -p "$myoutdir"
    
    configManta.py --bam={input.normal} --tumorBam={input.tumor} --referenceFasta {params.genome} --exome --runDir "$myoutdir"
    
    $myoutdir/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g 12
    """


rule annotsv:
    input: 
        sv_vcf=os.path.join(output_somatic_sv, "manta_out","{samples}","results/variants/somaticSV.vcf.gz")
    output: 
        annot_tsv = os.path.join(output_somatic_sv, "manta_out", "annotated_results","{samples}.annotsv.tsv"),
    params: 
        normalsample = lambda w: [pairs_dict[w.samples]],
        tumorsample = "{samples}",
        genome = config['references']['MAF_GENOME'],
        plot_script = config['scripts']['sv_plot'],
        rname = 'annotsv'
    threads: 4
    envmodules:
        'annotsv/2.2',
        'R/4.1'
    # container: config['images']['wes_base']
    shell: """
    myoutdir="$(dirname {output.annot_tsv})"
    if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi

    module load annotsv/2.2;
    AnnotSV -SVinputFile {input.sv_vcf} \\
            -genomeBuild {params.genome} \\
            -outputFile {output.annot_tsv} \\
            -outputDir $(dirname {output.annot_tsv})
    
    """
