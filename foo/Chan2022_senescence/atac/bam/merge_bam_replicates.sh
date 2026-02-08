
##############
## Settings ##
##############

INPUT_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac"
OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac/merged"

mkdir -p $OUTDIR; cd $OUTDIR

####################
## Samtools merge ##
####################

nohup samtools merge ${OUTDIR}/p8_Control.bam ${INPUT_DIR}/pA_Control_REP1.bam ${INPUT_DIR}/pA_Control_REP2.bam ${INPUT_DIR}/pA_Control_REP3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/p8_TG.bam ${INPUT_DIR}/pA_TG_REP1.bam ${INPUT_DIR}/pA_TG_REP2.bam ${INPUT_DIR}/pA_TG_REP3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/p11_Control.bam ${INPUT_DIR}/pB_Control_REP1.bam ${INPUT_DIR}/pB_Control_REP2.bam ${INPUT_DIR}/pB_Control_REP3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/p11_TG.bam ${INPUT_DIR}/pB_TG_REP1.bam ${INPUT_DIR}/pB_TG_REP2.bam ${INPUT_DIR}/pB_TG_REP3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/p14_Control.bam ${INPUT_DIR}/pC_Control_REP1.bam ${INPUT_DIR}/pC_Control_REP2.bam ${INPUT_DIR}/pC_Control_REP3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/p14_TG.bam ${INPUT_DIR}/pC_TG_REP1.bam ${INPUT_DIR}/pC_TG_REP2.bam ${INPUT_DIR}/pC_TG_REP3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/p19_Control.bam ${INPUT_DIR}/pD_Control_REP1.bam ${INPUT_DIR}/pD_Control_REP2.bam ${INPUT_DIR}/pD_Control_REP3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/p19_TG.bam ${INPUT_DIR}/pD_TG_REP1.bam ${INPUT_DIR}/pD_TG_REP2.bam ${INPUT_DIR}/pD_TG_REP3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/p22_Control.bam ${INPUT_DIR}/pE_Control_REP1.bam ${INPUT_DIR}/pE_Control_REP2.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/p22_TG.bam ${INPUT_DIR}/pE_TG_REP1.bam ${INPUT_DIR}/pE_TG_REP2.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/p24_Control.bam ${INPUT_DIR}/pF_Control_REP1.bam ${INPUT_DIR}/pF_Control_REP2.bam ${INPUT_DIR}/pF_Control_REP3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/p24_TG.bam ${INPUT_DIR}/pF_TG_REP1.bam ${INPUT_DIR}/pF_TG_REP2.bam ${INPUT_DIR}/pF_TG_REP3.bam -@ 1 & 

nohup samtools merge ${OUTDIR}/p29_Control.bam ${INPUT_DIR}/pG_Control_REP1.bam ${INPUT_DIR}/pG_Control_REP2.bam ${INPUT_DIR}/pG_Control_REP3.bam -@ 1 & 
nohup samtools merge ${OUTDIR}/p29_TG.bam ${INPUT_DIR}/pG_TG_REP1.bam ${INPUT_DIR}/pG_TG_REP2.bam -@ 1 & 