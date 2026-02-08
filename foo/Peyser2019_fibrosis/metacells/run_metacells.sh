conda activate metacells

SCRIPT="/share/rargelaguet/altos_github/fibrosis/Peyser2019_fibrosis/metacells/run_metacell.py"
ANNDATA="/share/rargelaguet/data/Peyser2019_fibrosis/processed/anndata.h5ad"
METADATA="/share/rargelaguet/data/Peyser2019_fibrosis/processed/cell_metadata.txt.gz"
OUTDIR="/share/rargelaguet/data/Peyser2019_fibrosis/results/metacells/per_condition"

mkdir -p $OUTDIR; cd $OUTDIR

samples=( "Saline Control" "Bleomycin Fibrosis Model" "Bleomycin Fibrosis plus Nintedanib" )

# i="FSFA"
for i in "${samples[@]}"; do
    out_i=$(echo "$i" | sed 's/ /-/g')
    if ! [ -f "${OUTDIR}/${i}/cell2metacell_assignment.txt.gz" ]; then
        echo $i
        nohup python $SCRIPT --anndata $ANNDATA --metadata $METADATA --outdir ${OUTDIR}/${out_i} --samples $i --percent_metacells 0.05 &> ${out_i}_nohup.out &
        # sbatch -n 4 --mem 15G -o slurm_${i}.out --wrap "python $SCRIPT --anndata $ANNDATA --metadata $METADATA --outdir ${OUTDIR}/${i} --samples $i --percent_metacells 0.05"
    fi
done
