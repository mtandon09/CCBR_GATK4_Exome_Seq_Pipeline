
rule sobdetect_get:
    input: 
    output: SOBDetector_jar=SOBDetector_JARFILE
    params: rname="get_sobdetector"
    shell:
      """
      wget https://github.com/mikdio/SOBDetector/releases/download/v1.0.2/SOBDetector_v1.0.2.jar -O {output.SOBDetector_jar};
      """
    
    
rule sobdetect_pass1:
    # input: vcf=os.path.join(output_somatic_snpindels,"{vc_outdir}","vcf","{samples}.FINAL.vcf"),
    input: vcf=os.path.join(output_somatic_snpindels,"{vc_outdir}","vcf","{samples}.FINAL.norm.vcf"),
           bam=os.path.join(output_bamdir,"final_bams","{samples}.bam"),
           SOBDetector_jar=SOBDetector_JARFILE
    output: pass1_vcf=os.path.join(SOBDetector_out,"{vc_outdir}","pass1","{samples}.sobdetect.vcf"),
            pass1_info=os.path.join(SOBDetector_out,"{vc_outdir}","pass1","{samples}.info")
    params: chrom=chroms, rname="sobdetect1"
    shell:"""
        module load samtools;
        myoutdir="$(dirname {output.pass1_vcf})"
        if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
        
        # touch $(echo "{input.bam}" | sed -e 's/.bam$/.bai/') ## SOBdetector sometimes complains about index being older than bam... Removing for now...
        
        java -jar {input.SOBDetector_jar} --input-type VCF --input-variants "{input.vcf}" --input-bam {input.bam} --output-variants {output.pass1_vcf} --only-passed false
        
        bcftools query -f '%INFO/numF1R2Alt\\t%INFO/numF2R1Alt\\t%INFO/numF1R2Ref\\t%INFO/numF2R1Ref\\t%INFO/numF1R2Other\\t%INFO/numF2R1Other\\t%INFO/SOB\\n' {output.pass1_vcf}| \
          awk '{{if ($1 != "."){{tum_alt=$1+$2; tum_depth=$1+$2+$3+$4+$5+$6; if (tum_depth==0){{tum_af=1}} else {{tum_af=tum_alt/tum_depth }}; print tum_alt,tum_depth,tum_af,$7}}}}' > {output.pass1_info} 
        """

    
rule sobdetect_cohort_params:
  input: info_files=expand(os.path.join(SOBDetector_out,"{{vc_outdir}}","pass1","{samples}.info"), samples=ffpe_sample_list)
  output: all_info_file=os.path.join(SOBDetector_out,"{vc_outdir}","pass1","all_samples.info"),
          params_file=os.path.join(SOBDetector_out,"{vc_outdir}","cohort_params.txt")
  params: rname="sobdetect_params"
  shell:
    """
    echo -e "#TUMOR.alt\\tTUMOR.depth\\tTUMOR.AF\\tSOB\\tFS\\tSOR\\tTLOD\\tReadPosRankSum" > {output.all_info_file}
    cat {input.info_files} >> {output.all_info_file}
    
    ## This is really stupid; I'm calculating mean and standard deviation in bash; a python one liner might be better?
    grep -v ^# {output.all_info_file} | \
      awk '{{ total1 += $1; ss1 += $1^2; total2 += $2; ss2 += $2^2; total3 += $3; ss3 += $3^2; total4 += $4; ss4 += $4^2 }} END {{ print total1/NR,total2/NR,total3/NR,total4/NR; print sqrt(ss1/NR-(total1/NR)^2),sqrt(ss2/NR-(total2/NR)^2),sqrt(ss3/NR-(total3/NR)^3),sqrt(ss4/NR-(total4/NR)^2) }}' > {output.params_file}
    """

    
rule sobdetect_pass2:
  # input: vcf=os.path.join(output_somatic_snpindels,"{vc_outdir}","vcf","{samples}.FINAL.vcf"),
  input: vcf=os.path.join(output_somatic_snpindels,"{vc_outdir}","vcf","{samples}.FINAL.norm.vcf"),
         bam=os.path.join(output_bamdir,"final_bams","{samples}.bam"),
         SOBDetector_jar=SOBDetector_JARFILE,
         params_file=os.path.join(SOBDetector_out,"{vc_outdir}","cohort_params.txt")
  output: pass2_vcf=os.path.join(SOBDetector_out,"{vc_outdir}","pass2","{samples}.sobdetect.vcf"),
          pass2_info=os.path.join(SOBDetector_out,"{vc_outdir}","pass2","{samples}.info"),
          filtered_vcf=os.path.join(SOBDetector_out,"{vc_outdir}","pass2","{samples}.artifact_filtered.vcf.gz")
  params: chrom=chroms,
          ver_bcftools=config['tools']['bcftools']['version'],
          rname="sobdetect2"
  shell:
    """
    module load samtools;
    date
    echo "Running SOBDetector..."
    myoutdir="$(dirname {output.pass2_vcf})"
    if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
    
    java -jar {input.SOBDetector_jar} --input-type VCF --input-variants "{input.vcf}" --input-bam "{input.bam}" --output-variants "{output.pass2_vcf}" --only-passed true --standardization-parameters "{input.params_file}"
    
    echo "Done with SOBDetector."
    date
    
    echo "Making info table..."
    module load bcftools/{params.ver_bcftools}
    bcftools query -f '%INFO/numF1R2Alt\\t%INFO/numF2R1Alt\\t%INFO/numF1R2Ref\\t%INFO/numF2R1Ref\\t%INFO/numF1R2Other\\t%INFO/numF2R1Other\\t%INFO/SOB\\n' "{output.pass2_vcf}" | \
      awk '{{if ($1 != "."){{tum_alt=$1+$2; tum_depth=$1+$2+$3+$4+$5+$6; if (tum_depth==0){{tum_af=1}} else {{tum_af=tum_alt/tum_depth }}; print tum_alt,tum_depth,tum_af,$7}}}}' > "{output.pass2_info}"
    echo "Done making info table.!"
    date

    echo "Filtering out artifacts..."
    
    if [ "{wildcards.vc_outdir}" == "{config[output_params][MERGED_SOMATIC_OUTDIR]}" ]; then
        echo "... and adding 'set' annotation back from merged variants..."
        bgzip --threads $((SLURM_CPUS_PER_TASK - 1)) -c "{input.vcf}" > "{input.vcf}.gz"
        bcftools index -f -t "{input.vcf}.gz"
        
        bgzip --threads $((SLURM_CPUS_PER_TASK - 1)) -c "{output.pass2_vcf}" > "{output.pass2_vcf}.gz"
        bcftools index -f -t "{output.pass2_vcf}.gz"
        
        bcftools annotate -a "{input.vcf}.gz" -c "INFO/set" -e 'INFO/pArtifact < 0.05' -Oz -o {output.filtered_vcf} {output.pass2_vcf}.gz
    else
        bcftools filter -e 'INFO/pArtifact < 0.05' -Oz -o {output.filtered_vcf} {output.pass2_vcf}
        bcftools index -f -t {output.filtered_vcf}
    fi
    
    echo "Done."
    date
    """

rule sobdetect_metrics:
  input: pass1_vcf=expand(os.path.join(SOBDetector_out,"{{vc_outdir}}","pass1","{samples}.sobdetect.vcf"), samples=ffpe_sample_list),
         pass2_vcf=expand(os.path.join(SOBDetector_out,"{{vc_outdir}}","pass2","{samples}.sobdetect.vcf"), samples=ffpe_sample_list)
  output: count_table=os.path.join(SOBDetector_out,"{vc_outdir}","metrics","variant_count_table.txt"),
          full_metric_table=os.path.join(SOBDetector_out,"{vc_outdir}","metrics","all_metrics.txt")
  params: ver_bcftools=config['tools']['bcftools']['version'],
          rname="sobdetect_metrics"
  shell:
    """
    echo -e "#ID\\tDefaultParam\\tCohortParam\\tTotalVariants" > {output.count_table}
    echo -e "#SAMPLE_ID\\tParam\\tCHROM\\tPOS\\tnumF1R2Alt\\tnumF2R1Alt\\tnumF1R2Ref\\tnumF2R1Ref\\tnumF1R2Other\\tnumF2R1Other\\tSOB\\tpArtifact\\tFS\\tSOR\\tTLOD\\tReadPosRankSum" > {output.full_metric_table}
    
    P1FILES=({input.pass1_vcf})
    P2FILES=({input.pass2_vcf})
    module load bcftools/{params.ver_bcftools}
    date
    for (( i=0; i<${{#P1FILES[@]}}; i++ )); do
      MYID=$(basename -s ".sobdetect.vcf" ${{P1FILES[$i]}})
      echo "Collecting metrics from $MYID..."
      total_count=$(grep -v ^# ${{P1FILES[$i]}} | wc -l)
      count_1p=$(bcftools query -f '%INFO/pArtifact\n' ${{P1FILES[$i]}} | awk '{{if ($1 != "." && $1 < 0.05){{print}}}}' | wc -l)
      count_2p=$(bcftools query -f '%INFO/pArtifact\n' ${{P2FILES[$i]}} | awk '{{if ($1 != "." && $1 < 0.05){{print}}}}' | wc -l)
      
      echo -e "$MYID\\t$count_1p\\t$count_2p\\t$total_count" >> {output.count_table}
      
      bcftools query -f '%CHROM\\t%POS\\t%INFO/numF1R2Alt\\t%INFO/numF2R1Alt\\t%INFO/numF1R2Ref\\t%INFO/numF2R1Ref\\t%INFO/numF1R2Other\\t%INFO/numF2R1Other\\t%INFO/SOB\\t%INFO/pArtifact\n' ${{P1FILES[$i]}} | awk -v id=$MYID 'BEGIN{{OFS="\t"}}{{print id,"PASS_1",$0}}' >> {output.full_metric_table}
      bcftools query -f '%CHROM\\t%POS\\t%INFO/numF1R2Alt\\t%INFO/numF2R1Alt\\t%INFO/numF1R2Ref\\t%INFO/numF2R1Ref\\t%INFO/numF1R2Other\\t%INFO/numF2R1Other\\t%INFO/SOB\\t%INFO/pArtifact\n' ${{P2FILES[$i]}} | awk -v id=$MYID 'BEGIN{{OFS="\t"}}{{print id,"PASS_2",$0}}' >> {output.full_metric_table}
    done
    
    echo "Done"
    date
    """

rule ffpefilter_mafs:
  input: filtered_vcf=os.path.join(SOBDetector_out,"{vc_outdir}","pass2","{samples}.artifact_filtered.vcf.gz"),
         vcf2maf_script=VCF2MAF_WRAPPER
  output: maf=os.path.join(output_somatic_base,SOBDetector_out,"{vc_outdir}","maf","{samples}.maf")
  params: tumorsample="{samples}",genome=config['references']['MAF_GENOME'],filtervcf=config['references']['MAF_FILTERVCF'],rname="pl:vcf2maf"
  shell:
    """
    date
    echo "Converting to MAF..."
    bash {input.vcf2maf_script} --vcf {input.filtered_vcf} --maf {output.maf} --tid {params.tumorsample} --genome {params.genome} --threads "$((SLURM_CPUS_PER_TASK-1))" --info "set"
    echo "Done converting to MAF..."
    date
    """
  
rule collect_ffpefilter_mafs:
  input: mafs=expand(os.path.join(SOBDetector_out,"{{vc_outdir}}","maf","{samples}"+".maf"), samples=ffpe_sample_list)
  output: maf=os.path.join(output_somatic_base,SOBDetector_out,"{vc_outdir}","cohort_summary","all_somatic_variants.maf")
  params: rname="combine_maf"
  shell:"""    
    date
    echo "Combining MAFs..."
    awk 'FNR==1 && NR!=1 {{ while (/^#/||/^Hugo_Symbol/) getline; }} 1 {{print}}' {input.mafs} > {output.maf}
    echo "Done..."
    date
  """