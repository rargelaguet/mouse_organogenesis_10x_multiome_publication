snakemake --cores 1
snakemake --cores 1 --dry-run -p
snakemake --cores 10 -j 99 --latency-wait 90 -p --cluster "sbatch -n {threads} --mem {resources.mem_mb}M"

