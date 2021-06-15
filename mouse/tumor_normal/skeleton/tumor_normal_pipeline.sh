
#!/usr/bin/bash
module load snakemake/5.24.1; snakemake --stats snakemake.stats --restart-times 0 --rerun-incomplete --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e slurmfiles/slurm-%j_{params.rname}.out -o slurmfiles/slurm-%j_{params.rname}.out" --cluster-config cluster.json --keep-going --snakefile tumor_normal_hg38.snakemake  --config "input_params={'BAM_SOURCE':'/data/tandonm/pl_test_data/human/bams','PAIRS_FILE':'pairs.tsv','BASE_OUTDIR':'/scratch/tandonm/pipeline_out/tn_test_1'}" -j 500
 > snakemake.log 2>&1
