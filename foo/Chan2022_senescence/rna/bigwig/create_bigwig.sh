
##############
## Settings ##
##############

BAMCOVERAGE="/share/rargelaguet/mambaforge/envs/ucsc/bin/bamCoverage"
BAM_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/rna/merged"
OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bigwig/rna/merged"

mkdir -p $OUTDIR; cd $OUTDIR

# samples=( "pA_control_rep1" "pA_control_rep2" "pA_control_rep3" "pA_tg_rep1" "pA_tg_rep2" "pA_tg_rep3" "pB_control_rep1" "pB_control_rep2" "pB_control_rep3" "pB_tg_rep1" "pB_tg_rep2" "pB_tg_rep3" "pC_control_rep1" "pC_control_rep2" "pC_control_rep3" "pC_tg_rep1" "pC_tg_rep2" "pC_tg_rep3" "pD_control_rep1" "pD_control_rep2" "pD_control_rep3" "pD_tg_rep1" "pD_tg_rep2" "pD_tg_rep3" "pE_control_rep1" "pE_control_rep2" "pE_control_rep3" "pE_tg_rep1" "pE_tg_rep2" "pE_tg_rep3" "pF_control_rep1" "pF_control_rep2" "pF_control_rep3" "pF_tg_rep1" "pF_tg_rep2" "pF_tg_rep3" "pG_control_rep1" "pG_control_rep2" "pG_control_rep3" "pG_tg_rep1" "pG_tg_rep2" "pG_tg_rep3" )
samples=( "Control_p11" "Control_p14" "Control_p19" "Control_p22" "Control_p24" "Control_p29" "Control_p8" "TG_p11" "TG_p14" "TG_p19" "TG_p22" "TG_p24" "TG_p29" "TG_p8" )

cores=1

################################
## create index for bam files ##
################################

cd $BAM_DIR

for i in "${samples[@]}"; do
	if ! [ -f "${i}.bam.bai" ]; then
		echo $i
	    nohup samtools index ${i}.bam & 
	fi
done

####################
## create bigwigs ##
####################

mkdir -p $OUTDIR; cd $OUTDIR

for i in "${samples[@]}"; do
	if ! [ -f "${OUTDIR}/${i}.bw" ]; then
		echo $i
		nohup $BAMCOVERAGE -b ${BAM_DIR}/${i}.bam -o $OUTDIR/$i.bw --binSize 25 --numberOfProcessors $cores --normalizeUsing CPM &
    fi
done