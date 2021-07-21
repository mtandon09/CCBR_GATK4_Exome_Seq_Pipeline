
rule trimmomatic:
    input:  r1=os.path.join(input_fqdir, "{samples}.R1.fastq.gz"),
            r2=os.path.join(input_fqdir, "{samples}.R2.fastq.gz")
    output: one=temp(os.path.join(output_fqdir,"{samples}.R1.trimmed.fastq.gz")),
            two=temp(os.path.join(output_fqdir,"{samples}.R1.trimmed.unpair.fastq.gz")),
            three=temp(os.path.join(output_fqdir,"{samples}.R2.trimmed.fastq.gz")),
            four=temp(os.path.join(output_fqdir,"{samples}.R2.trimmed.unpair.fastq.gz")),
            err=os.path.join(output_fqdir,"{samples}_run_trimmomatic.err")
    params: adapterfile=config['references']['trimmomatic.adapters'],ver=config['tools']['trimmomatic']['version'],rname="pl:trimmomatic"
    # threads: 32
    shell:  """
            module load trimmomatic/{params.ver};
            myoutdir="$(dirname {output.one})"
            if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
            
            trimmomatic PE -threads $((SLURM_CPUS_PER_TASK-1)) -phred33 {input.r1} {input.r2} {output.one} {output.two} {output.three} {output.four} ILLUMINACLIP:{params.adapterfile}:3:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:20 2> {output.err}
            
            """

rule bwa_mem:
    input:  os.path.join(output_fqdir,"{samples}.R1.trimmed.fastq.gz"),
            os.path.join(output_fqdir,"{samples}.R2.trimmed.fastq.gz")
    output: temp(os.path.join(output_bamdir,"preprocessing","{samples}.raw_map.bam"))
    params: genome=config['references']['GENOME'],sample = "{samples}",ver_samtools=config['tools']['samtools']['version'],ver_bwa=config['tools']['bwa']['version'],rname="pl:bwamem"
    # threads: 32
    shell:  """
            module load samtools/{params.ver_samtools}
            module load bwa/{params.ver_bwa}
            myoutdir="$(dirname {output})"
            if [ ! -d "$myoutdir" ]; then mkdir -p "$myoutdir"; fi
            bwa mem -M -R \'@RG\\tID:{params.sample}\\tSM:{params.sample}\\tPL:illumina\\tLB:{params.sample}\\tPU:{params.sample}\\tCN:hgsc\\tDS:wes\' -t $((SLURM_CPUS_PER_TASK-1)) {params.genome} {input} | /usr/local/apps/samblaster/0.1.25/bin/samblaster -M | samtools sort -@12 -m 4G - -o {output}
            """

rule raw_index:
      input:  bam=os.path.join(output_bamdir,"preprocessing","{samples}.raw_map.bam")
      output: bai=temp(os.path.join(output_bamdir,"preprocessing","{samples}.raw_map.bai")),
      params: ver_samtools=config['tools']['samtools']['version'],rname="raw_index"
      shell:  """
              module load samtools/{params.ver_samtools}
              samtools index -@ 2 {input.bam} {output.bai}
              """
rule realign:
      input:  bam=os.path.join(output_bamdir,"preprocessing","{samples}.raw_map.bam"),
              bai=os.path.join(output_bamdir,"preprocessing","{samples}.raw_map.bai"),
      output: bam=temp(os.path.join(output_bamdir,"preprocessing","{samples}.realign.bam")),
              int=temp(os.path.join(output_bamdir,"preprocessing","{samples}.intervals")),
      params: genome=config['references']['GENOME'],knowns=config['references']['KNOWNINDELS'],ver_gatk=config['tools']['gatk3']['version'],rname="realign"
      shell:  """
              module load GATK/{params.ver_gatk}
              java -Xmx48g -jar $GATK_JAR -T RealignerTargetCreator -I {input.bam} -R {params.genome} -o {output.int} {params.knowns}
              java -Xmx48g -jar $GATK_JAR -T IndelRealigner -R {params.genome} -I {input.bam} {params.knowns} --use_jdk_inflater --use_jdk_deflater -targetIntervals {output.int} -o {output.bam}
              """

rule gatk_recal:
      input:  os.path.join(output_bamdir,"preprocessing","{samples}.realign.bam")
      output: bam=os.path.join(input_bamdir,"{samples}.input.bam"),
              re=temp(os.path.join(output_bamdir,"preprocessing","{samples}_recal_data.grp"))
      params: genome=config['references']['GENOME'],knowns=config['references']['KNOWNRECAL'],ver_gatk=config['tools']['gatk4']['version'],rname="recal"
      threads: 24
      shell:  """
              module load GATK/{params.ver_gatk}
              gatk --java-options '-Xmx48g' BaseRecalibrator --input {input} --reference {params.genome} {params.knowns} --output {output.re}
              gatk --java-options '-Xmx48g' ApplyBQSR --reference {params.genome} --input {input} --bqsr-recal-file {output.re} --output {output.bam} --use-jdk-inflater --use-jdk-deflater
              """

rule bam_check:
      input:  bam=os.path.join(input_bamdir,"{samples}.input.bam")
      output: bam=os.path.join(output_bamdir,"final_bams","{samples}.bam"),
              bai=os.path.join(output_bamdir,"final_bams","{samples}.bai"),
              bai2=os.path.join(output_bamdir,"final_bams","{samples}.bam.bai"),
      params: ver_samtools=config['tools']['samtools']['version'],
              ver_gatk=config['tools']['gatk4']['version'],
              rname="bam_check"
      shell:  """
              module load samtools/{params.ver_samtools}
              module load GATK/{params.ver_gatk}
              
              sample={wildcards.samples}
              ID=$sample
              PL="ILLUMINA"  ### Should be exposed as a config param
              LB="na"        ### Should be exposed as a config param (?)
              
              ##Check if there is no header or any of the info
              HEADER=`samtools view -H {input.bam} | grep ^@RG`
              if [[ "$HEADER" != "" ]]; then
                  t=(${{HEADER//\t/ }})
                  echo ${{t[1]}}
                  ID=`printf '%s\n' "${{t[@]}}" | grep -P '^ID' | cut -d":" -f2` #(${{t[1]//:/ }})
                  PL=`printf '%s\n' "${{t[@]}}" | grep -P '^PL' | cut -d":" -f2` #(${{t[3]//:/ }})
                  LB=`printf '%s\n' "${{t[@]}}" | grep -P '^LB' | cut -d":" -f2` #(${{t[2]//:/ }})
                  if [[ "$ID" == "$sample" ]]; then
                      echo "The header of the BAM file is correct"
                  else
                      ID=$sample
                  fi
              fi
              gatk AddOrReplaceReadGroups --INPUT {input.bam} --OUTPUT {output.bam} --RGID ${{ID}} --RGLB ${{LB}} --RGPL ${{PL}} --RGSM ${{ID}} --RGPU na
              
              samtools index -@ 2 {output.bam} {output.bai}
              cp {output.bai} {output.bai2}
              """
