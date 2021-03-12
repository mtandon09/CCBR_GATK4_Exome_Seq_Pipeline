#!/usr/bin/bash

pipeline_dir="../skeleton/"

pngdir="$(pwd)/rulegraphs"  ## Where to output rule graph PNGs
if [ ! -d $pngdir ]; then mkdir -p $pngdir; fi

cd $pipeline_dir

## Default (uses params defined in the references_hg38.json)
./run.sh --rulegraph "$pngdir/rules.default.png"

## Start from BAM files
# Since the config file already contains a source fastq parameter, need to override it with a non-existent filepath
./run.sh --rulegraph "$pngdir/rules.fromBAM.png" --sourcefq "foo" --sourcebam "/data/tandonm/pl_test_data/human/bams"

## Add CNV calling
# Should be set of one of 'true', 't', or 'yes' (case-insensitive)
./run.sh --rulegraph "$pngdir/rules.CNV.png" --sourcefq "foo" --cnv "True"

## Add FFPE filtering with SOBDetector
# Should be set of one of 'true', 't', or 'yes' (case-insensitive)
./run.sh --rulegraph "$pngdir/rules.FFPE.png" --sourcefq "foo" --ffpe "True"

## Custom set of variant callers
#  This doesn't work yet; needs to be deleted from the config json!!
custom_config_str="output_params={'SOMATIC_VCF':{'mutect2':'','mutect':''}}"
./run.sh --rulegraph "$pngdir/rules.custom_callers.png" --config "$custom_config_str"

