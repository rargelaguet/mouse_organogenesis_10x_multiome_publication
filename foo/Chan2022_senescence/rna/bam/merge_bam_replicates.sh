
##############
## Settings ##
##############

INPUT_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/rna"
OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/rna/merged"

mkdir -p $OUTDIR; cd $OUTDIR

####################
## Samtools merge ##
####################

nohup samtools merge ${OUTDIR}/Control_p8.bam ${INPUT_DIR}/pA_control_rep1.bam ${INPUT_DIR}/pA_control_rep2.bam ${INPUT_DIR}/pA_control_rep3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/TG_p8.bam ${INPUT_DIR}/pA_tg_rep1.bam ${INPUT_DIR}/pA_tg_rep2.bam ${INPUT_DIR}/pA_tg_rep3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/Control_p11.bam ${INPUT_DIR}/pB_control_rep1.bam ${INPUT_DIR}/pB_control_rep2.bam ${INPUT_DIR}/pB_control_rep3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/TG_p11.bam ${INPUT_DIR}/pB_tg_rep1.bam ${INPUT_DIR}/pB_tg_rep2.bam ${INPUT_DIR}/pB_tg_rep3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/Control_p14.bam ${INPUT_DIR}/pC_control_rep1.bam ${INPUT_DIR}/pC_control_rep2.bam ${INPUT_DIR}/pC_control_rep3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/TG_p14.bam ${INPUT_DIR}/pC_tg_rep1.bam ${INPUT_DIR}/pC_tg_rep2.bam ${INPUT_DIR}/pC_tg_rep3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/Control_p19.bam ${INPUT_DIR}/pD_control_rep1.bam ${INPUT_DIR}/pD_control_rep2.bam ${INPUT_DIR}/pD_control_rep3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/TG_p19.bam ${INPUT_DIR}/pD_tg_rep1.bam ${INPUT_DIR}/pD_tg_rep2.bam ${INPUT_DIR}/pD_tg_rep3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/Control_p22.bam ${INPUT_DIR}/pE_control_rep1.bam ${INPUT_DIR}/pE_control_rep2.bam ${INPUT_DIR}/pE_control_rep3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/TG_p22.bam ${INPUT_DIR}/pE_tg_rep1.bam ${INPUT_DIR}/pE_tg_rep2.bam ${INPUT_DIR}/pE_tg_rep3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/Control_p24.bam ${INPUT_DIR}/pF_control_rep1.bam ${INPUT_DIR}/pF_control_rep2.bam ${INPUT_DIR}/pF_control_rep3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/TG_p24.bam ${INPUT_DIR}/pF_tg_rep1.bam ${INPUT_DIR}/pF_tg_rep2.bam ${INPUT_DIR}/pF_tg_rep3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/Control_p29.bam ${INPUT_DIR}/pG_control_rep1.bam ${INPUT_DIR}/pG_control_rep2.bam ${INPUT_DIR}/pG_control_rep3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/TG_p29.bam ${INPUT_DIR}/pG_tg_rep1.bam ${INPUT_DIR}/pG_tg_rep2.bam ${INPUT_DIR}/pG_tg_rep3.bam -@ 1 & 
