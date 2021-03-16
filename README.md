# hg38 pipeline
-------------------------------------------------
~~Currently only the **tumor_normal** pipeline is up-to-date. I'm still toying with the idea of making a unified pipeline that will run tumor-only if no pairs info is provided.~~
Both `tumor_only` and `tumor_normal` are functional.  It got too complicated to merge them together, so they are separate for now.

### Deploying
So far, I've usually been creating a copy of the skeleton for each run. Then I either edit the config json (`references_hg38.json` *prob wanna change this filename* lol) and call snakemake manually. Or increasingly, I've built up the `run.sh` script to handle most of the things I needed during development.  Theoretically by setting the input/output options correctly, a common skeleton can be used for any number of jobs but I have not tested it thoroughly.

The defaults to watch out for or set explicitly when calling `run.sh`:
- `pairs.tsv` is read from the skeleton folder (`--pairs` arugment)
- Slurm output files are stored in the skeleton (`--slurmfiles` argument)
- The snakemake log file for each submission is stored in `submit.log` in the skeleton (currently not configurable)

### Config tests
Here's some bash commmands for deploying this on Biowulf and running the config tests described [here](https://github.com/mtandon09/exome_pipeline_dev_mt/tree/main/tumor_normal/config_tests).
```
## Download the repo
git clone https://github.com/mtandon09/exome_pipeline_dev_mt.git

## Set up a working directory
WORKDIR="/scratch/$USER/wes_pipe_test"
if [ ! -d $WORKDIR ]; then mkdir -p "$WORKDIR"; fi

## Copy over the skeleton
cp -r exome_pipeline_dev_mt/tumor_normal/* $WORKDIR

## Try to run the config tests
cd $WORKDIR/config_tests

./make_rulegraphs.sh

```

### `run.sh`
This is meant to be a catch-all driver script to do common tasks for running the pipeline.

```
usage: run.sh [-h] [--sourcefq SOURCEFQ] [--sourcebam SOURCEBAM]
              [--pairs PAIRS] [--callers CALLERS] [--targets TARGETS]
              [--ffpe FFPE] [--cnv CNV] [--outdir OUTDIR] [--dryrun DRYRUN]
              [--unlock UNLOCK] [--until UNTIL] [--local LOCAL]
              [--slurmdir SLURMDIR] [--rulegraph RULEGRAPH] [--config CONFIG]

Run muh pipelinezz

optional arguments:
  -h, --help            show this help message and exit
  --sourcefq SOURCEFQ   [input_params] Path to directory containing paired
                        FASTQ files
  --sourcebam SOURCEBAM
                        [input_params] Path to directory containing paired BAM
                        files. If '--sourcefq' is also defined, the sample IDs
                        should match the FASTQ files.
  --pairs PAIRS         [input_params] TSV file containing two columns with
                        tumor and normal sample IDs, one pair per line. The
                        header needs to be 'Tumor' for the tumor column and
                        'Normal' for the normal column.
  --callers CALLERS     [input_params] list of mutation callers, comma-
                        separated in single quotes. Default:
                        "'mutect2','mutect','strelka','vardict','varscan'"
  --targets TARGETS     [input_params] Path to exome targets BED file
  --ffpe FFPE           [input_params] Add FFPE filtering step (set to one of
                        'true', 't', or 'yes' (case-insensitive))
  --cnv CNV             [input_params] Add CNV calling step (set to one of
                        'true', 't', or 'yes' (case-insensitive))
  --outdir OUTDIR       [input_params] Location to store pipeline output
  --dryrun DRYRUN       Dry-run only (provide any non-empty string)
  --unlock UNLOCK       Unlock working directory (provide any non-empty
                        string)
  --until UNTIL         Rule name to stop at; passed to snakemake's '--until'
                        argument
  --local LOCAL         Number of jobs to run in parallel locally; does not
                        submit to slurm, so only use on an interactive node
  --slurmdir SLURMDIR   Path to output slurm files to
  --rulegraph RULEGRAPH
                        Path to a PNG file to which the rules DAG will be
                        written
  --config CONFIG       Manually set the 'input_params' section of the
                        snakemake config file. Overrides any [input_params]
                        arguments.
```

## Example Run
The test data contains 10 million randomly sampled read pairs from three human cell lines.

Starting from BAM files, the basic pipeline (germline and somatic SNP calling) for this test data runs in a little over 2 hours total. 

[Tumor-Normal Snakemake Report](tumor_normal/skeleton/report.html)

```

## Run this in the skeleton directory (or a copy)
test_data_fq="/data/tandonm/pl_test_data/human/fastq"
test_data_bam="/data/tandonm/pl_test_data/human/bams"
pairs_file="pairs.tsv"

## See if the dry-run works
./run.sh --dryrun 1 \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam"

## Without the --dryrun argument, the job will be submitted
./run.sh --pairs "$pairs_file" \
         --sourcebam "$test_data_bam"

```


