#!/usr/bin/bash
module load snakemake; snakemake --stats snakemake.stats --restart-times 0 --rerun-incomplete --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e slurmfiles/slurm-%j_{params.rname}.out -o slurmfiles/slurm-%j_{params.rname}.out" --cluster-config cluster.json --keep-going --snakefile tumor_only_hg38.snakemake  -j 500
 > snakemake.log 2>&1
