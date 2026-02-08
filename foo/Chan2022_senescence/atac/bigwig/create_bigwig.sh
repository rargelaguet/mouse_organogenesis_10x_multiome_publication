conda activate ucsc

##############
## Settings ##
##############

BAM_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac/merged"
OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/bigwig"
cores=1

samples=( "Control_p11" "Control_p14" "Control_p19" "Control_p22" "Control_p24" "Control_p29" "Control_p8" "TG_p11" "TG_p14" "TG_p19" "TG_p22" "TG_p24" "TG_p29" "TG_p8" )

################################
## create index for bam files ##
################################

cd $BAM_DIR

for i in "${samples[@]}"; do
	if ! [ -f "${i}.bam.bai" ]; then
	    nohup samtools index ${i}.bam &
	fi
done

####################
## create bigwigs ##
####################

mkdir -p $OUTDIR; cd $OUTDIR

# i="Control_p11"
for i in "${samples[@]}"; do
	if ! [ -f "${OUTDIR}/${i}.bw" ]; then
		echo $i
		nohup bamCoverage -b ${BAM_DIR}/${i}.bam -o $OUTDIR/$i.bw --binSize 25 --numberOfProcessors $cores --normalizeUsing CPM &
    fi
done

