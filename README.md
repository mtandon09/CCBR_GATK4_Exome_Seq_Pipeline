# hg38 WES pipeline

The `unified` folder contains the latest workflow containing all the changes in the separate `tumor_normal` and `tumor_only` workflows from before.

Essentially, all you need to start the pipeline is a folder of gzip-ed fastq files (and a target bed file if not using Agilent SureSelect Human All Exon v7). This will run germline and somatic calling in tumor-only mode.  If a pairs file is provided, it will be run in tumor-normal mode.

## Example Run
The test data contains 10 million randomly sampled read pairs from three human cell lines.

Starting from FASTQ files, the basic pipeline (germline and somatic SNP calling) for this test data runs in about 4 hours. 

Here's some bash commmands for deploying this on Biowulf.

```
## Download the repo
git clone https://github.com/mtandon09/exome_pipeline_dev_mt.git
cd exome_pipeline_dev_mt/unified/skeleton

## Define some parameters
output_dir="/scratch/$USER/wes_pipe_test"
test_data_fq="/data/tandonm/pl_test_data/human/fastq"

## See if the dry-run works
./run.sh --sourcefq "$test_data_fq" \
         --outdir "$output_dir" \
         --dryrun 1

## Without the --dryrun argument, the job will be submitted
./run.sh --sourcefq "$test_data_fq" --outdir "$output_dir"

```


## Customizing run parameters
### `run.sh`
This is meant to be a catch-all driver script to do common tasks for running the pipeline.

```
usage: run.sh [-h] [--sourcefq SOURCEFQ] [--sourcebam SOURCEBAM] [--mode MODE]
              [--pairs PAIRS] [--callers CALLERS] [--targets TARGETS]
              [--ffpe FFPE] [--cnv CNV] [--outdir OUTDIR] [--dryrun DRYRUN]
              [--unlock UNLOCK] [--until UNTIL] [--local LOCAL]
              [--slurmdir SLURMDIR] [--rulegraph RULEGRAPH] [--report REPORT]
              [--config CONFIG]

Tumor-only Pipeline

optional arguments:
  -h, --help            show this help message and exit
  --sourcefq SOURCEFQ   [input_params] Path to directory containing paired
                        FASTQ files
  --sourcebam SOURCEBAM
                        [input_params] Path to directory containing paired BAM
                        files. If '--sourcefq' is also defined, the sample IDs
                        should match the FASTQ files.
  --mode MODE           [input_params] One of 'tumor_only' or 'paired'
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
  --report REPORT       Path to an HTML file to which the snakemake report
                        will be written
  --config CONFIG       Manually set the 'input_params' section of the
                        snakemake config file. Overrides any [input_params]
                        arguments.
```


### Notes/To-do

- Currently this works best if you copy the skeleton for each run. Mostly because the snakemake log is written to the skeleton (I guess wherever `run.sh` is called from)
- Need a bam to fastq rule; currently, running from bam will cause a couple of the QC rules that require fastqs to fail
- Replace run.sh with a python program; ideally generate config.json on-the-fly 

