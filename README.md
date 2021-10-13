# CCBR Exome-Seq Pipeline for Biowulf

# Usage

### Clone the repo
```
git clone https://github.com/mtandon09/exome_pipeline_dev_mt.git
```

### Set up a new run
```
## The folder with the pipeline skeleton
PIPE_DIR="$(pwd)/exome_pipeline_dev_mt/pipeline"

## New folder where we want to run it
NEW_FOLDER="/data/$USER/exome_test"

## Create the new folder
if [ ! -d $NEW_FOLDER ]; then mkdir -p $NEW_FOLDER; fi

## Create a subdirectory to store the pipeline in the new folder
cd $NEW_FOLDER
mkdir pipeline

## Copy over the pipeline skeleton
cp -r $PIPE_DIR/* pipeline/

cd pipeline

## Can run this to see all options for the run command
./run.sh --help


## Set up variables for input params for a run
# Folder of fastq.gz files
fq_dir="/data/tandonm/pl_test_data/human/fastq"

# Location of subdirectory to hold pipeline output
output_dir="$NEW_FOLDER/output"

# Location of targets bed file (hg38)
# This is the default file used in the config file, so no need to specify unless it's different
bed_file="/data/CCBR_Pipeliner/db/PipeDB/lib/Agilent_SSv7_allExons_hg38.bed"

# Pairs file; must contain column header "Tumor" and "Normal" (in any order) listing sample IDs
pairs_file="/data/tandonm/pl_test_data/human/pairs"


## Try a dry run to see if it all rules are compiled successfully
./run.sh --sourcefq "$fq_dir" --outdir "$output_dir" --targets "$bed_file" --pairs "$pairs_file" --ffpe "True" --dryrun 1

## Submit to the cluster by excluding the --dryrun flag
./run.sh --sourcefq "$fq_dir" --outdir "$output_dir" --targets "$bed_file" --pairs "$pairs_file" --ffpe "True"

```




