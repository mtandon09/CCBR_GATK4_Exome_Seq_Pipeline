#!/usr/bin/bash


pipeline_dir="../skeleton/"
test_data_fq="/data/tandonm/pl_test_data/human/fastq"
test_data_bam="/data/tandonm/pl_test_data/human/bams"
pairs_file="pairs.tsv"

pngdir="$(pwd)/rulegraphs"  ## Where to output rule graph PNGs
if [ ! -d $pngdir ]; then mkdir -p $pngdir; fi

cd $pipeline_dir

## Start from Fastq
./run.sh --rulegraph "$pngdir/rules.fromFQ.png" \
         --pairs "$pairs_file" \
         --sourcefq "$test_data_fq"

## Start from BAM files
./run.sh --rulegraph "$pngdir/rules.fromBAM.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam"

## Add CNV calling
# Should be set of one of 'true', 't', or 'yes' (case-insensitive)
./run.sh --rulegraph "$pngdir/rules.CNV.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --cnv "True"

## Add FFPE filtering with SOBDetector
# Should be set of one of 'true', 't', or 'yes' (case-insensitive)
./run.sh --rulegraph "$pngdir/rules.FFPE.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --ffpe "True"

## Custom set of variant callers
#  Each caller should be in single-quotes with commas between callers
#  Available options: 'mutect2','mutect','strelka','vardict','varscan'
caller_str="'strelka','varscan'"
./run.sh --rulegraph "$pngdir/rules.custom_callers.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --callers "$caller_str"

# Note that if a single caller is selected, the merge_somatic_callers step will be automatically omitted
./run.sh --rulegraph "$pngdir/rules.single_caller.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --callers "'vardict'"

# This command will try a dry-run, only selecting mutect2, and only until merge_somatic_callers
./run.sh --dryrun 1 \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --callers "'mutect2'" \
         --until merge_somatic_callers

# Snakemake should find that no jobs need to be run
#[+] Loading snakemake  5.24.1 
#Building DAG of jobs...
#Nothing to be done.