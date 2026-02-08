conda activate metacells

SCRIPT="/share/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/metacells/run/run_metacell.py"
ANNDATA="/share/rargelaguet/data/Xing2020_bj_oskm/processed/rna/scTE/anndata.h5ad"
METADATA="/share/rargelaguet/data/Xing2020_bj_oskm/processed/rna/scTE/cell_metadata.txt.gz"
OUTDIR="/share/rargelaguet/data/Xing2020_bj_oskm/results/rna/metacells/scTE"

mkdir -p $OUTDIR; cd $OUTDIR

samples=( "sample-SRR7657442_D0_10X_merged" "sample-SRR7657443_D2_10X_merged" "sample-SRR7657444_D8_10X_merged" "sample-SRR7657445_D12_10X_merged" "sample-SRR7657446_D16__10X_merged" "sample-SRR7657447_D16__10X_merged" )

# i="sample-SRR7657442_D0_10X_merged"
for i in "${samples[@]}"; do
    if ! [ -f "${OUTDIR}/${i}/cell2metacell_assignment.txt.gz" ]; then
        # python run_metacell.py --anndata $ANNDATA --metadata $METADATA --outdir ${OUTDIR}/${i} --samples $i --percent_metacells 0.05
        nohup python $SCRIPT --anndata $ANNDATA --metadata $METADATA --outdir ${OUTDIR}/${i} --samples $i --percent_metacells 0.05 &> ${i}_nohup.out &
        # sbatch -n 4 --mem 30G -o slurm_${i}.out --wrap "python $SCRIPT --anndata $ANNDATA --metadata $METADATA --outdir ${OUTDIR}/${i} --samples $i --percent_metacells 0.05"
    fi
done
