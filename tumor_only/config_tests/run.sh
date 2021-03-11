#!/usr/bin/bash

source /data/CCBR_Pipeliner/db/PipeDB/lib/vcf2maf_resources/scripts/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--dryrun',required=False, default='', help='Dry-run only')
parser.add_argument('--unlock',required=False, default='', help='Unlock working directory (pass any non-empty string to this option)')
parser.add_argument('--until',required=False,  default='', help='Rule name to stop at; passed to snakemake\'s \'--until\' argument')
parser.add_argument('--local',required=False,  default='', help='Number of jobs to run in parallel locally; does not submit to slurm, so need to be on an interactive node')
parser.add_argument('--slurmdir',required=False, default='slurmfiles', help='Path to output slurm files to')
parser.add_argument('--rulegraph',required=False, default='', help='Path to a PNG file to which the rules DAG will be written')
EOF


SNAKEFILE="tumor_only_hg38.snakemake"
if [ ! -z $UNLOCK ]; then
    module load snakemake
    snakemake -j2 --unlock --snakefile $SNAKEFILE --rerun-incomplete
    exit 0
fi

if [ ! -z $RULEGRAPH ]; then
    module load snakemake
    module load graphviz
    snakemake --snakefile $SNAKEFILE --rerun-incomplete --rulegraph | dot -Tpng > $RULEGRAPH
    echo "Wrote PNG to $RULEGRAPH."
    exit 0
fi

untilarg=""
if [ ! -z $UNTIL ]; then untilarg="--until $UNTIL"; fi

if [ ! -d $SLURMDIR ]; then mkdir -p $SLURMDIR; fi

CLUSTER_OPTS="sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e $SLURMDIR/slurm-%j_{params.rname}.out -o $SLURMDIR/slurm-%j_{params.rname}.out";

dryrunarg=""
if [ ! -z $DRYRUN ]; then dryrunarg="-n"; fi

#smk_cmd_base="snakemake --stats snakemake.stats --restart-times 0 --rerun-incomplete -j 500 --cluster "$CLUSTER_OPTS" --cluster-config cluster.json --keep-going --snakefile $SNAKEFILE $untilarg $dryrunarg > snakemake.log 2>&1;"
smk_cmd_base="module load snakemake; snakemake --stats snakemake.stats --restart-times 0 --rerun-incomplete --cluster \"$CLUSTER_OPTS\" --cluster-config cluster.json --keep-going --snakefile $SNAKEFILE $untilarg"


NJOBS="500"
if [ ! -z $DRYRUN ]; then
    eval "$smk_cmd_base -npr -j $NJOBS"
else
    if [ ! -z $LOCAL ]; then
        echo "module load snakemake; snakemake --rerun-incomplete --snakefile $SNAKEFILE -j $LOCAL $untilarg"
        eval "module load snakemake; snakemake --rerun-incomplete --snakefile $SNAKEFILE -j $LOCAL $untilarg"
    else
        echo -e "#!/usr/bin/bash\n$smk_cmd_base -j $NJOBS\n > snakemake.log 2>&1" > tumor_only_pipeline.sh
        echo "Submitting pipeline to cluster... "
        sbatch --cpus-per-task=2 --mem=12g --time 5-00:00:00 --partition ccr,norm --output submit.log --error submit.log tumor_only_pipeline.sh
    fi
fi
