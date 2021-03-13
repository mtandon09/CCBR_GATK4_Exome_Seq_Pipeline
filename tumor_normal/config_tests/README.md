# Configurable pipeline steps

The configurable parameters are part of the `input_params` section of the [config JSON](https://github.com/mtandon09/exome_pipeline_dev_mt/blob/main/config_tests/references_hg38.json).
```
...
"input_params": {
      "FASTQ_SOURCE": "/data/tandonm/pl_test_data/human/fastq",
      "BAM_SOURCE": "/data/tandonm/pl_test_data/human/bams",
      "PAIRS_FILE": "pairs.tsv",
      "EXOME_TARGETS": "/data/CCBR_Pipeliner/db/PipeDB/lib/Agilent_SSv7_allExons_hg38.bed",
      "FFPE_FILTER": "false",
      "CNV_CALLING": "false",
      "BASE_OUTDIR": "tn_out_1"
    },
...
```
We can also control these parameters with the `run.sh` script in the pipeline skeleton.  You can run `./run.sh --help` to display all available options.

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

At a minimum, these arguments must be provided
- `--pairs`
- `--targets`
- At least one of `--sourcefq` or `--sourcebam`


FASTQs for testing can be found here: `/data/tandonm/pl_test_data/human/fastq`
BAMs for testing can be found here: `/data/tandonm/pl_test_data/human/bams`

```
pipeline_dir="../skeleton/"
test_data_fq="/data/tandonm/pl_test_data/human/fastq"
test_data_bam="/data/tandonm/pl_test_data/human/bams"
pairs_file="pairs.tsv"

pngdir="$(pwd)/rulegraphs"  ## Where to output rule graph PNGs
```

---------------------------------------------------
## Start from FASTQ files

This rulegraph can be generated with the following command on Biowulf:
```
## Start from Fastq
./run.sh --rulegraph "$pngdir/rules.fromFQ.png" \
         --pairs "$pairs_file" \
         --sourcefq "$test_data_fq"
```

![Start from fastq files](rulegraphs/rules.fromFQ.png)



---------------------------------------------------
## Start from BAM files


If the fastq directory does not contain fastq files, the BAM directory will be used automatically.
```
## Start from BAM files
./run.sh --rulegraph "$pngdir/rules.fromBAM.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam"
```

![Start from BAM files](rulegraphs/rules.fromBAM.png)




---------------------------------------------------
## Add copy number variant (CNV) calling

This is turned off by default, so use the option in `run.sh` to turn it on.
```
## Add CNV calling
## Should be set of one of 'true', 't', or 'yes' (case-insensitive)
./run.sh --rulegraph "$pngdir/rules.CNV.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --cnv "True"
```

![With CNV calling](rulegraphs/rules.CNV.png)


---------------------------------------------------
## Add FFPE artifact filtering

Currently using [`SOBDetector`](https://github.com/mikdio/SOBDetector) to flag FFPE artifacts for somatic calls only. 

This is turned off by default, so use the option in `run.sh` to turn it on.
```
## Add FFPE filtering with SOBDetector
## Should be set of one of 'true', 't', or 'yes' (case-insensitive)
./run.sh --rulegraph "$pngdir/rules.FFPE.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --ffpe "True"
```

![With SOBDetector](rulegraphs/rules.FFPE.png)


---------------------------------------------------
## Turn variant callers on/off

The default config json file will run all six variant valling steps
- `mutect2`
- `mutect` (v 1.XX)
- `strelka`
- `vardict`
- `varscan`
- Merged calls from all callers

You can customize this with the `--callers` command-line option.  Each caller should be in single-quotes with commas between callers.
```
## Custom set of variant callers
##  Available options: 'mutect2','mutect','strelka','vardict','varscan'
caller_str="'strelka','varscan'"
./run.sh --rulegraph "$pngdir/rules.custom_callers.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --callers "$caller_str"
```

![Custom caller set](rulegraphs/rules.custom_callers.png)


Note that if a single caller is selected, the `merge_somatic_callers` step will be automatically omitted

```
./run.sh --rulegraph "$pngdir/rules.single_caller.png" \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --callers "'vardict'"
```

![Single caller](rulegraphs/rules.single_caller.png)

And if you try to dry-run to that rule like this:
```
./run.sh --dryrun 1 \
         --pairs "$pairs_file" \
         --sourcebam "$test_data_bam" \
         --callers "'mutect2'" \
         --until merge_somatic_callers
```

Snakemake will determine that no jobs should be run.
```
[+] Loading snakemake  5.24.1 
Building DAG of jobs...
Nothing to be done.
```


