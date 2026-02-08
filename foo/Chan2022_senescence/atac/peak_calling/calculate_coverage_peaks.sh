conda activate ucsc

##############
## Settings ##
##############

BED_FILE="/share/rargelaguet/data/BJ_passaging_CHAMP/results/atac/peak_calling/macs3/merged_peaks.bed"

BAM_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac"
samples=( "pA_Control_REP1" "pA_Control_REP2" "pA_Control_REP3" "pA_TG_REP1" "pA_TG_REP2" "pA_TG_REP3" "pB_Control_REP1" "pB_Control_REP2" "pB_Control_REP3" "pB_TG_REP1" "pB_TG_REP2" "pB_TG_REP3" "pC_Control_REP1" "pC_Control_REP2" "pC_Control_REP3" "pC_TG_REP1" "pC_TG_REP2" "pC_TG_REP3" "pD_Control_REP1" "pD_Control_REP2" "pD_Control_REP3" "pD_TG_REP1" "pD_TG_REP2" "pD_TG_REP3" "pE_Control_REP1" "pE_Control_REP2" "pE_TG_REP1" "pE_TG_REP2" "pF_Control_REP1" "pF_Control_REP2" "pF_Control_REP3" "pF_TG_REP1" "pF_TG_REP2" "pF_TG_REP3" "pG_Control_REP1" "pG_Control_REP2" "pG_Control_REP3" "pG_TG_REP1" "pG_TG_REP2" )

# BAM_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac/merged"
# samples=( "Control_p11" "Control_p14" "Control_p19" "Control_p22" "Control_p24" "Control_p29" "Control_p8" "TG_p11" "TG_p14" "TG_p19" "TG_p22" "TG_p24" "TG_p29" "TG_p8" )

OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/coverage_peaks"

mkdir -p $OUTDIR; cd $OUTDIR

#######################
## bedtools coverage ##
#######################

# # i="Control_p11"
# for i in "${samples[@]}"; do
#     echo $i
#     nohup bedtools coverage -d -a ${BED_FILE} -b ${BAM_DIR}/${i}.bam > ${OUTDIR}/${i}.bed &
#     # sbatch -n 4 --mem 15G -o slurm_${i}.out --wrap "bedtools coverage -d -a ${BED_FILE} -b ${BAM_DIR}/${i}.bam"
# done

nohup bedtools multicov -bams ${BAM_DIR}/*.bam -bed $BED_FILE > ${OUTDIR}/atac_coverage_peaks.count &