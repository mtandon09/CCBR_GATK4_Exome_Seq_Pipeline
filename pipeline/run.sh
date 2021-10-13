#!/usr/bin/bash
set -e

source /data/CCBR_Pipeliner/db/PipeDB/lib/vcf2maf_resources/scripts/argparse.bash || exit 1
ARGPARSE_DESCRIPTION="Tumor-only Pipeline"
argparse "$@" <<EOF || exit 1
parser.add_argument('--sourcefq',required=False, default='', help='[input_params] Path to directory containing paired FASTQ files')
parser.add_argument('--sourcebam',required=False, default='', help='[input_params] Path to directory containing paired BAM files.  If \'--sourcefq\' is also defined, the sample IDs should match the FASTQ files.')
parser.add_argument('--mode',required=False, default='tumor_only', help='[input_params] One of \'tumor_only\' or \'paired\'')
parser.add_argument('--pairs',required=False, default='', help='[input_params] TSV file containing two columns with tumor and normal sample IDs, one pair per line.  The header needs to be \'Tumor\' for the tumor column and \'Normal\' for the normal column.')
parser.add_argument('--callers',required=False, default='', help='[input_params] list of mutation callers, comma-separated in single quotes. Default: "\'mutect2\',\'mutect\',\'strelka\',\'vardict\',\'varscan\'"')
parser.add_argument('--targets',required=False, default='', help='[input_params] Path to exome targets BED file')
parser.add_argument('--ffpe',required=False, default='', help='[input_params] Add FFPE filtering step (set to one of \'true\', \'t\', or \'yes\' (case-insensitive))')
parser.add_argument('--cnv',required=False, default='', help='[input_params] Add CNV calling step (set to one of \'true\', \'t\', or \'yes\' (case-insensitive))')
parser.add_argument('--outdir',required=False, default='pipeline_output', help='[input_params] Location to store pipeline output')
parser.add_argument('--dryrun',required=False, default='', help='Dry-run only (provide any non-empty string)')
parser.add_argument('--unlock',required=False, default='', help='Unlock working directory (provide any non-empty string)')
parser.add_argument('--until',required=False,  default='', help='Rule name to stop at; passed to snakemake\'s \'--until\' argument')
parser.add_argument('--local',required=False,  default='', help='Number of jobs to run in parallel locally; does not submit to slurm, so only use on an interactive node')
parser.add_argument('--slurmdir',required=False, default='', help='Path to output slurm files to')
parser.add_argument('--rulegraph',required=False, default='', help='Path to a PNG file to which the rules DAG will be written')
parser.add_argument('--report',required=False, default='', help='Path to an HTML file to which the snakemake report will be written')
parser.add_argument('--config',required=False, default='', help='Manually set the \'input_params\' section of the snakemake config file. Overrides any [input_params] arguments.')
EOF

SNAKEFILE="workflow/Snakefile"
CLUSTER_JSON="config/cluster.json"
SUBMIT_SCRIPT="submit_pipe.sh"
SMK_VERSION="6.5.3"

untilarg=""
if [ ! -z $UNTIL ]; then untilarg="--until $UNTIL"; fi

MODE="'TN_MODE':'$MODE'"

configarg=""
if [ ! -z $CONFIG ];
    then configarg="--config \"$CONFIG\"";
else
    if [ ! -z $SOURCEFQ ]; then SOURCEFQ="'FASTQ_SOURCE':'$SOURCEFQ'"; fi
    if [ ! -z $SOURCEBAM ]; then SOURCEBAM="'BAM_SOURCE':'$SOURCEBAM'"; fi
    if [ ! -z $PAIRS ]; then PAIRS="'PAIRS_FILE':'$PAIRS'"; MODE="'TN_MODE':'paired'"; fi
    if [ ! -z $CALLERS ]; then CALLERS="'VARIANT_CALLERS':[$CALLERS]"; fi
    if [ ! -z $TARGETS ]; then TARGETS="'EXOME_TARGETS':'$TARGETS'"; fi
    if [ ! -z $FFPE ]; then FFPE="'FFPE_FILTER':'$FFPE'"; fi
    if [ ! -z $CNV ]; then CNV="'CNV_CALLING':'$CNV'"; fi
    if [ ! -z $OUTDIR ]; then CNV="'BASE_OUTDIR':'$OUTDIR'"; fi
    inputparams=($SOURCEFQ $SOURCEBAM $MODE $PAIRS $CALLERS $TARGETS $FFPE $CNV)
    for i in ${inputparams[@]}; do configarg="$configarg,$i"; done
    if [ ! -z $configarg ]; then configarg="--config \"input_params={$(echo -e $configarg | sed -e 's/^,//')}\""; fi
fi
if [ ! -z $UNLOCK ]; then
    module load snakemake/$SMK_VERSION
    echo $configarg
    mycmd="snakemake -j2 --snakefile $SNAKEFILE --rerun-incomplete $configarg --unlock"
    eval $mycmd
    exit 0
fi
if [ ! -z $RULEGRAPH ]; then
    module load snakemake/$SMK_VERSION
    module load graphviz
    mycmd="snakemake --snakefile $SNAKEFILE --rerun-incomplete --rulegraph $configarg | dot -Tpng > $RULEGRAPH"
    #echo $mycmd
    eval $mycmd
    echo "Wrote rule graph to $RULEGRAPH. Exiting..."
    exit 0
fi
if [ ! -z $REPORT ]; then
    module load snakemake/$SMK_VERSION
    module load graphviz
    mycmd="snakemake --snakefile $SNAKEFILE --rerun-incomplete --report $REPORT $configarg"
    #echo $mycmd
    eval $mycmd
    echo "Wrote report to $REPORT. Exiting..."
    exit 0
fi

if [ -z $SLURMDIR ]; then SLURMDIR="$OUTDIR/slurmfiles"; fi
if [ ! -d $SLURMDIR ]; then mkdir -p $SLURMDIR; fi

CLUSTER_OPTS="sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e $SLURMDIR/slurm-%j_{params.rname}.out -o $SLURMDIR/slurm-%j_{params.rname}.out";

dryrunarg=""
if [ ! -z $DRYRUN ]; then dryrunarg="-n"; fi

#smk_cmd_base="snakemake --stats snakemake.stats --restart-times 0 --rerun-incomplete -j 500 --cluster "$CLUSTER_OPTS" --cluster-config cluster.json --keep-going --snakefile $SNAKEFILE $untilarg $dryrunarg > snakemake.log 2>&1;"
smk_cmd_base="module load snakemake/$SMK_VERSION; snakemake --stats snakemake.stats --rerun-incomplete --cluster \"$CLUSTER_OPTS\" --cluster-config $CLUSTER_JSON --keep-going --snakefile $SNAKEFILE $untilarg $configarg"

## Can get group names with this, could use it to build a "--groups" or "--group-components" argument
# grep -r "group: " workflow/* | sed -r 's/.*group: "(.*)"/\1/' | sort | uniq

NJOBS="500"
if [ ! -z $DRYRUN ]; then
    #echo "$smk_cmd_base -npr -j $NJOBS"
    eval "$smk_cmd_base -npr -j $NJOBS"
else
    if [ ! -z $LOCAL ]; then
        echo "module load snakemake/$SMK_VERSION; module load graphviz; snakemake --rerun-incomplete --snakefile $SNAKEFILE -j $LOCAL $untilarg $configarg"
        eval "module load snakemake/$SMK_VERSION; module load graphviz; snakemake --rerun-incomplete --snakefile $SNAKEFILE -j $LOCAL $untilarg $configarg"
    else
        echo -e "#!/usr/bin/bash\n$smk_cmd_base -j $NJOBS  --local-cores 12 --restart-times 1 --latency-wait 120 \n &> snakemake.log" > $SUBMIT_SCRIPT
        echo "Submitting pipeline to cluster... "
        primaryID=$(sbatch --cpus-per-task=32 --mem=96g --time 5-00:00:00 --gres=lscratch:500 --partition ccr,norm --output submit.log --error submit.log $SUBMIT_SCRIPT)
        echo "Primary Job ID: $primaryID"
        
    fi
fi
