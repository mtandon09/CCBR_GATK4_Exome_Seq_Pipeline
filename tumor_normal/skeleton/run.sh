#!/usr/bin/bash
set -e

source /data/CCBR_Pipeliner/db/PipeDB/lib/vcf2maf_resources/scripts/argparse.bash || exit 1
ARGPARSE_DESCRIPTION="Run muh pipelinezz"
argparse "$@" <<EOF || exit 1
parser.add_argument('--sourcefq',required=False, default='', help='[input_params] Path to directory containing paired FASTQ files')
parser.add_argument('--sourcebam',required=False, default='', help='[input_params] Path to directory containing paired BAM files.  If \'--sourcefq\' is also defined, the sample IDs should match the FASTQ files.')
parser.add_argument('--pairs',required=False, default='', help='[input_params] TSV file containing two columns with tumor and normal sample IDs, one pair per line.  The header needs to be \'Tumor\' for the tumor column and \'Normal\' for the normal column.')
parser.add_argument('--targets',required=False, default='', help='[input_params] Path to exome targets BED file')
parser.add_argument('--ffpe',required=False, default='', help='[input_params] Add FFPE filtering step (set to one of \'true\', \'t\', or \'yes\' (case-insensitive))')
parser.add_argument('--cnv',required=False, default='', help='[input_params] Add CNV calling step (set to one of \'true\', \'t\', or \'yes\' (case-insensitive))')
parser.add_argument('--outdir',required=False, default='', help='[input_params] Location to store pipeline output')
parser.add_argument('--dryrun',required=False, default='', help='Dry-run only (provide any non-empty string)')
parser.add_argument('--unlock',required=False, default='', help='Unlock working directory (provide any non-empty string)')
parser.add_argument('--until',required=False,  default='', help='Rule name to stop at; passed to snakemake\'s \'--until\' argument')
parser.add_argument('--local',required=False,  default='', help='Number of jobs to run in parallel locally; does not submit to slurm, so only use on an interactive node')
parser.add_argument('--slurmdir',required=False, default='slurmfiles', help='Path to output slurm files to')
parser.add_argument('--rulegraph',required=False, default='', help='Path to a PNG file to which the rules DAG will be written')
parser.add_argument('--config',required=False, default='', help='Manually set the \'input_params\' section of the snakemake config file. Overrides any [input_params] arguments.')
EOF

SNAKEFILE="tumor_normal_hg38.snakemake"

if [ ! -z $UNLOCK ]; then
    module load snakemake/5.24.1
    snakemake -j2 --unlock --snakefile $SNAKEFILE --rerun-incomplete
    exit 0
fi

untilarg=""
if [ ! -z $UNTIL ]; then untilarg="--until $UNTIL"; fi


configarg=""
if [ ! -z $CONFIG ];
    then configarg="--config \"$CONFIG\"";
else
    if [ ! -z $SOURCEFQ ]; then SOURCEFQ="'FASTQ_SOURCE':'$SOURCEFQ'"; fi
    if [ ! -z $SOURCEBAM ]; then SOURCEBAM="'BAM_SOURCE':'$SOURCEBAM'"; fi
    if [ ! -z $PAIRS ]; then PAIRS="'PAIRS_FILE':'$PAIRS'"; fi
    if [ ! -z $TARGETS ]; then TARGETS="'EXOME_TARGETS':'$TARGETS'"; fi
    if [ ! -z $FFPE ]; then FFPE="'FFPE_FILTER':'$FFPE'"; fi
    if [ ! -z $CNV ]; then CNV="'CNV_CALLING':'$CNV'"; fi
    if [ ! -z $OUTDIR ]; then CNV="'BASE_OUTDIR':'$OUTDIR'"; fi
    inputparams=($SOURCEFQ $SOURCEBAM $PAIRS $TARGETS $FFPE $CNV)
    for i in ${inputparams[@]}; do configarg="$configarg,$i"; done
    if [ ! -z $configarg ]; then configarg="--config \"input_params={$(echo -e $configarg | sed -e 's/^,//')}\""; fi
fi

if [ ! -z $RULEGRAPH ]; then
    module load snakemake/5.24.1
    module load graphviz
    mycmd="snakemake --snakefile $SNAKEFILE --rerun-incomplete --rulegraph $configarg | dot -Tpng > $RULEGRAPH"
    echo $mycmd
    eval $mycmd
    echo "Wrote rule graph to $RULEGRAPH. Exiting..."
    exit 0
fi

if [ ! -d $SLURMDIR ]; then mkdir -p $SLURMDIR; fi

CLUSTER_OPTS="sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e $SLURMDIR/slurm-%j_{params.rname}.out -o $SLURMDIR/slurm-%j_{params.rname}.out";

dryrunarg=""
if [ ! -z $DRYRUN ]; then dryrunarg="-n"; fi

#smk_cmd_base="snakemake --stats snakemake.stats --restart-times 0 --rerun-incomplete -j 500 --cluster "$CLUSTER_OPTS" --cluster-config cluster.json --keep-going --snakefile $SNAKEFILE $untilarg $dryrunarg > snakemake.log 2>&1;"
smk_cmd_base="module load snakemake/5.24.1; snakemake --stats snakemake.stats --restart-times 0 --rerun-incomplete --cluster \"$CLUSTER_OPTS\" --cluster-config cluster.json --keep-going --snakefile $SNAKEFILE $untilarg $configarg"


NJOBS="500"
if [ ! -z $DRYRUN ]; then
    echo "$smk_cmd_base -npr -j $NJOBS"
    eval "$smk_cmd_base -npr -j $NJOBS"
else
    if [ ! -z $LOCAL ]; then
        echo "module load snakemake/5.24.1; snakemake --rerun-incomplete --snakefile $SNAKEFILE -j $LOCAL $untilarg"
        eval "module load snakemake/5.24.1; snakemake --rerun-incomplete --snakefile $SNAKEFILE -j $LOCAL $untilarg"
    else
        echo -e "#!/usr/bin/bash\n$smk_cmd_base -j $NJOBS\n > snakemake.log 2>&1" > tumor_normal_pipeline.sh
        echo "Submitting pipeline to cluster... "
        sbatch --cpus-per-task=2 --mem=12g --time 5-00:00:00 --partition ccr,norm --output submit.log --error submit.log tumor_normal_pipeline.sh
    fi
fi
