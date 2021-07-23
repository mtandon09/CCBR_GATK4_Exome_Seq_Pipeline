#!/usr/bin/bash

pipeline_dir="../skeleton/"
test_data_fq="/data/tandonm/pl_test_data/human/fastq"
#test_data_bam="/data/tandonm/pl_test_data/human/bams"

pngdir="$(pwd)/rulegraphs"  ## Where to output rule graph PNGs
if [ ! -d $pngdir ]; then mkdir -p $pngdir; fi

cd $pipeline_dir

## TUMOR ONLY (default), Start from Fastq
./run.sh --rulegraph "$pngdir/rules.fromFQ.png" \
         --sourcefq "$test_data_fq"

## Start from BAM files
#./run.sh --rulegraph "$pngdir/rules.fromBAM.png" \
#         --sourcebam "$test_data_bam"

## TUMOR ONLY (default), Add FFPE filtering with SOBDetector
# Should be set of one of 'true', 't', or 'yes' (case-insensitive)
./run.sh --rulegraph "$pngdir/rules.FFPE.png" \
         --sourcefq "$test_data_fq" \
         --ffpe "True"

## TUMOR ONLY (default), Custom set of variant callers
#  Each caller should be in single-quotes with commas between callers
#  Available options: 'mutect2','mutect','vardict','varscan'
caller_str="'mutect2','vardict'"
./run.sh --rulegraph "$pngdir/rules.custom_callers.png" \
         --sourcefq "$test_data_fq" \
         --callers "$caller_str"

# Note that if a single caller is selected, the merge_somatic_callers step will be automatically omitted
./run.sh --rulegraph "$pngdir/rules.single_caller.png" \
         --sourcefq "$test_data_fq" \
         --callers "'vardict'"

# This command will try a dry-run, only selecting mutect, and only until merge_somatic_callers
./run.sh --dryrun 1 \
         --sourcefq "$test_data_fq" \
         --callers "'mutect'" \
         --until "somatic_merge_callers"

# Snakemake should find that no jobs need to be run
# (because if only a single caller is selected, there is no need to merge)
#[+] Loading snakemake  5.24.1 
#Building DAG of jobs...
#Nothing to be done.



# A couple of illustrations of using paired mode
pairs_file="config/pairs.mixed.tsv"

## Literally everything possible
./run.sh --rulegraph "$pngdir/rules.paired_everything.png" \
         --sourcefq "$test_data_fq" \
         --pairs "$pairs_file" \
         --callers "'mutect2','strelka','vardict','mutect','varscan'" \
         --ffpe "True" \
         --cnv "True"

## Minimum possible?
./run.sh --rulegraph "$pngdir/rules.paired_minimum.png" \
         --sourcefq "$test_data_fq" \
         --pairs "$pairs_file" \
         --callers "'strelka'" \
         --until "multiqc"