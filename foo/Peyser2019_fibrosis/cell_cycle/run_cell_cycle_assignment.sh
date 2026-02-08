conda activate R_4.3

SCRIPT="/share/rargelaguet/altos_github/aging/fibros_reik_lab/rna/cell_cycle/cell_cycle_assignment_seurat.R"
OUTDIR="/share/rargelaguet/data/fibros_reik_lab/results/rna/cell_cycle"

samples=( "FSFA" "FSFB" "25F" "28F" "52F" "54F" "63F" "70M" )

# mkdir -p $OUTDIR; cd $OUTDIR
cd /share/rargelaguet/altos_github/aging/fibros_reik_lab/rna/cell_cycle

# i="FSFA"
for i in "${samples[@]}"; do
    if ! [ -f "${OUTDIR}/${i}_cell_cycle.txt.gz" ]; then
        echo $i
        nohup Rscript $SCRIPT --sample $i &> ${i}_nohup.out &
        # sbatch -n 4 --mem 15G -o slurm_${i}.out --wrap "Rscript $SCRIPT --sample $i"
    fi
done
