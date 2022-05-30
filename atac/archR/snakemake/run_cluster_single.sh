snakemake --cores 1 -j 1 --latency-wait 90 -p --cluster "sbatch -n {threads} --mem {resources.mem_mb}M"
