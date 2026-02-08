conda activate scanpy

SCRIPT="/share/rargelaguet/altos_github/reprogramming/Xing2020_bj_oskm/dimensionality_reduction/dimred_per_stage/dimred_scanpy_args.py"
ANNDATA="/share/rargelaguet/data/Xing2020_bj_oskm/processed/cellranger/anndata.h5ad"
METADATA="/share/rargelaguet/data/Xing2020_bj_oskm/processed/cellranger/cell_metadata.txt.gz"
OUTDIR="/share/rargelaguet/data/Xing2020_bj_oskm/results/dimensionality_reduction"

mkdir -p $OUTDIR; cd $OUTDIR

# samples=( "D00_Dox_1" "D00_Dox_2" "D00_Dox_Recover1" "D00_Dox_Recover2" "D00_noDox_1" "D00_noDox_2" "D02_Dox_1" "D02_Dox_2" "D02_Dox_Recover1" "D02_Dox_Recover2" "D04_Dox_1" "D04_Dox_2" "D04_Dox_Recover1" "D04_Dox_Recover2" "D04_noDox_1" "D04_noDox_2" "D06_Dox_1" "D06_Dox_2" "D06_Dox_Recover1" "D06_Dox_Recover2" "D06_Pulse_1" "D06_Pulse_2" "D06_Pulse_Recover1" "D06_Pulse_Recover2" "D08_Dox_1" "D08_Dox_2" "D08_Dox_Recover1" "D08_Dox_Recover2" "D08_noDox_1" "D08_noDox_2" "D12_noDox_1" "D12_noDox_2" )
groups=( "D00_noDox" "D04_noDox" "D08_noDox" "D12_noDox" "D00_Dox" "D02_Dox" "D04_Dox" "D06_Dox" "D08_Dox" "D00_Dox_Recover" "D02_Dox_Recover" "D04_Dox_Recover" "D06_Dox_Recover" "D08_Dox_Recover" "D06_Pulse" "D06_Pulse_Recover" )

# i="D00_noDox"
for i in "${groups[@]}"; do
    if ! [ -f "${OUTDIR}/${i}/X_pca.txt" ]; then
        # python $SCRIPT --anndata $ANNDATA --metadata $METADATA --outdir ${OUTDIR}/${i} --samples $i
        nohup python $SCRIPT --anndata $ANNDATA --metadata $METADATA --outdir ${OUTDIR}/${i} --samples $i &> ${i}_nohup.out &
        # sbatch -n 4 --mem 30G -o slurm_${i}.out --wrap "python $SCRIPT --anndata $ANNDATA --metadata $METADATA --outdir ${OUTDIR}/${i} --samples $i"
    fi
done