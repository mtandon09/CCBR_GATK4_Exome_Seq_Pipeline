#!/bin/bash
module load singularity
WORKDIR=
#. /opt2/conda/etc/profile.d/conda.sh
#conda activate python2

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="vcf2maf wrapper script"
#source /opt2/argparse.bash || exit 1
source /data/CCBR_Pipeliner/db/PipeDB/lib/vcf2maf_resources/scripts/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--vcf',required=True, help='input vcf file')
parser.add_argument('--maf',required=True, help='output maf file')
parser.add_argument('--genome',required=True, help='hg19/hg38/mm10')
parser.add_argument('--tid',required=False, help='TumorID in VCF ... sample id of tumor sample [first column name of VCF used by default]')
parser.add_argument('--nid',required=False, help='NormalID in VCF ... sample id of normal sample ... do not provide in case of tumor-only')
parser.add_argument('--maf_tid',required=False, default='',help='TumorID for MAF... output sample id of tumor sample [--tid used by defaul]')
parser.add_argument('--maf_nid',required=False, default='',help='NormalID for MAF ... output sample id of normal sample [--nid used by defaul]')
parser.add_argument('--info', required=False, default='', help='Comma-delimited names of INFO fields to retain as extra columns in MAF [default \'\']')
parser.add_argument('--addchr',required=False, default=True, help='Whether or not to add a chr prefix to chromosome names [default=True]')
parser.add_argument('--threads',required=False, default=0, help='Number of forks to use for VEP [default=1, or $SLURM_CPUS_PER_TASK-1 if on Biowulf]')
parser.add_argument('--vepresourcebundlepath',required=False, default='/data/CCBR_Pipeliner/db/PipeDB/lib/vcf2maf_resources',help='location of the resources required by vep')
EOF

ncbi_build="unknown"
species="unknown"
if [ $GENOME == "hg38" ]; then
        ncbi_build="GRCh38"
        species="homo_sapiens"
fi

if [ $GENOME == "hg19" ]; then
        ncbi_build="GRCh37"
        species="homo_sapiens"
fi

if [ $GENOME == "mm10" ]; then
        ncbi_build="GRCm38"
        species="mus_musculus"
fi

fa="/vepresourcebundlepath/fastas/${species}/${ncbi_build}.fa"
dotvep="/vepresourcebundlepath/VEP_tarballs/.vep"
filtervcf="/vepresourcebundlepath/filtervcf/${species}/${ncbi_build}.filter.vcf.gz"

OUTPUT_DIR=$(dirname $MAF)

## Detect file extension and unzip if necessary
VCF_BASENAME=$(basename $VCF)
vcf_file_extension=${VCF_BASENAME##*.}
if [ $vcf_file_extension == "gz" ]; then
        ungz_vcf="$OUTPUT_DIR/$(basename -s .gz $VCF)"
        echo "Decompressing VCF to $ungz_vcf..."
        gzip -dc $VCF > $ungz_vcf
        VCF=$ungz_vcf
fi
INPUT_DIR=$(dirname $VCF)

## This seems to solve this error from vcf2maf for certain VCFs:
## "ERROR: Your VCF uses CR line breaks, which we can't support. Please use LF or CRLF"
dos2unix -ascii $VCF

## Add chr prefix if requested and missing
chr_text="chr"
if [ $GENOME == "hg19" ]; then
        chr_text=""
else
        chr_text="chr"
fi
echo "Adding text '$chr_text' to chrom names..."
awk -v pre=$chr_text '{out_text=$0; gsub("^chr","",out_text); if( /^[^#]/ ) { out_text=pre out_text}; print out_text }' $VCF > $VCF.fixed
mv $VCF.fixed $VCF
#if [ ${ADDCHR:0:1}=="T" ]; then
#        echo "Adding chr prefix to chroms if missing..."
#        mv $VCF $VCF.original
#        awk '{out_text=$0; gsub("^chr","",out_text); if( /^[^#]/ ) { out_text="chr"out_text }; print out_text }' $VCF.original > $VCF
#fi

## Detect sample names in VCF
module load bcftools
VCF_SAMPLE_IDS=($(bcftools query -l $VCF))

## By default the first column will be the tumor (in case of single column)
TID_IDX=0
VCF_NID=""
NSAMPLES=${#VCF_SAMPLE_IDS[@]}
if [ $NSAMPLES -gt 1 ]; then
## just assume two columns, first normal then tumor
        NID_IDX=0
        TID_IDX=1 
        ## And assign tumor/normal IDs 
        ## Look through column names and see if they match provided IDs
        for (( i = 0; i < $NSAMPLES; i++ )); do
                echo "${VCF_SAMPLE_IDS[$i]}"
                if [ "${VCF_SAMPLE_IDS[$i]}" == "$TID" ]; then
                        TID_IDX=$i
                fi
                
                if [ "${VCF_SAMPLE_IDS[$i]}" == "$NID" ]; then
                        NID_IDX=$i
                fi
        done
        VCF_NID=${VCF_SAMPLE_IDS[$NID_IDX]}
fi
VCF_TID=${VCF_SAMPLE_IDS[$TID_IDX]}


if [ $MAF_TID == ""]; then
        MAF_TID=$VCF_TID
fi
if [ $MAF_NID == ""]; then
        MAF_NID=$VCF_NID
fi


if [ $THREADS==0 ]; then
        THREADS=1
        ## Use threads if on Biowulf
        if [ ! -z $SLURM_CPUS_PER_TASK ]; then
                THREADS=$((SLURM_CPUS_PER_TASK-1))
        fi
fi

#VEP_VERSION="97"
#module load VEP/$VEP_VERSION
#VEPRESOURCEBUNDLEPATH="${VEP_CACHEDIR}/${species}/${VEP_VERSION}_${ncbi_build}"

#OUTPUT_DIR=$(dirname $MAF)

PYTHONNOUSERSITE=1 singularity exec --cleanenv \
-B $INPUT_DIR:/data/vcfdir/,$VEPRESOURCEBUNDLEPATH/:/vepresourcebundlepath/,$OUTPUT_DIR:/mnt/outdir \
/data/CCBR_Pipeliner/db/PipeDB/db/SingularityImages/ccbr_vcf2maf_v0.1.0.sif \
vcf2maf.pl \
--vep-forks "$THREADS" \
--input-vcf "/data/vcfdir/$(basename $VCF)" \
--vcf-tumor-id "$VCF_TID" \
--vcf-normal-id "$VCF_NID" \
--output-maf "/mnt/outdir/$(basename $MAF)" \
--tumor-id "$MAF_TID" \
--normal-id "$MAF_NID" \
--vep-path "/opt/vep/src/ensembl-vep" \
--vep-data "$dotvep" \
--filter-vcf "$filtervcf" \
--ncbi-build "$ncbi_build" \
--species "$species" \
--retain-info "$INFO" \
--ref-fasta "$fa"

