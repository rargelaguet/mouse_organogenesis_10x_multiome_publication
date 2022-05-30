#################
## Run locally ##
#################

snakemake --use-conda --cores 1
snakemake --forceall --use-conda --cores 1
snakemake --forceall --use-conda --cores 1 --dry-run

#################################
## Run on the Babraham cluster ##
#################################

sbatch -n 1 --mem 5G snakemake --forceall -j 4 --use-conda --latency-wait 90 --cluster "sbatch -n 1 --mem 5G"
snakemake --forceall -j 4 --use-conda --latency-wait 90 --cluster "sbatch -n 1 --mem 12G"
snakemake -j 4 --use-conda --latency-wait 90 --cluster "sbatch -n 1 --mem {params.memory}G"