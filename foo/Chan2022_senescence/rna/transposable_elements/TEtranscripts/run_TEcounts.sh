conda activate TEtranscripts

##############
## Settings ##
##############

gtf_transposable_elements="/share/rargelaguet/data/hg38_regulation/transposable_elements/GRCh38_Ensembl_rmsk_TE.gtf.ind"
gtf_genes="/share/rargelaguet/data/ensembl/human/v109/GTF/Homo_sapiens.GRCh38.109.gtf.ind"

samples=( "p8_Control_REP1" "p8_Control_REP2" "p8_Control_REP3" "p8_TG_REP1" "p8_TG_REP2" "p8_TG_REP3" "p11_Control_REP1" "p11_Control_REP2" "p11_Control_REP3" "p11_TG_REP1" "p11_TG_REP2" "p11_TG_REP3" "p14_Control_REP1" "p14_Control_REP2" "p14_Control_REP3" "p14_TG_REP1" "p14_TG_REP2" "p14_TG_REP3" "p19_Control_REP1" "p19_Control_REP2" "p19_Control_REP3" "p19_TG_REP1" "p19_TG_REP2" "p19_TG_REP3" "p22_Control_REP1" "p22_Control_REP2" "p22_Control_REP3" "p22_TG_REP1" "p22_TG_REP2" "p22_TG_REP3" "p24_Control_REP1" "p24_Control_REP2" "p24_Control_REP3" "p24_TG_REP1" "p24_TG_REP2" "p24_TG_REP3" "p29_Control_REP1" "p29_Control_REP2" "p29_Control_REP3" "p29_TG_REP1" "p29_TG_REP2" "p29_TG_REP3" )

################################
## create index for GTF files ##
################################

# TEtranscripts_indexer=/share/rargelaguet/git/TEtranscripts/bin/TEtranscripts_indexer
# gtf_transposable_elements="/share/rargelaguet/data/hg38_regulation/transposable_elements/GRCh38_Ensembl_rmsk_TE.chr.gtf.ind"
# gtf_genes="/share/rargelaguet/data/ensembl/human/v109/GTF/Homo_sapiens.GRCh38.109.chr.gtf.ind"

# cd /share/rargelaguet/data/ensembl/human/v109/GTF/
# nohup $TEtranscripts_indexer --afile $gtf_genes --itype gene &

# cd /share/rargelaguet/data/hg38_regulation/transposable_elements/
# nohup $TEtranscripts_indexer --afile $gtf_transposable_elements --itype TE &

################################
## create index for bam files ##
################################

# cd /bi/group/reik/ricard/data/DellaValle2022_line1_progeria/star/first_trial

# for i in "${samples[@]}"; do
# 	if ! [ -f "${i}_Aligned.sortedByCoord.out.bam.bai" ]; then
# 		echo $i
# 	    job 10 1 "samtools index ${i}_Aligned.sortedByCoord.out.bam"
# 	fi
# done

#############
## TEcount ##
#############

BAMDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/rna/symbolic"
OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/rna/TEtranscripts"

mkdir -p $OUTDIR; cd $OUTDIR

# i="SRR23630177_GSM7067645_SKIN_AGE85_M_GR5_Gt125"
for i in "${samples[@]}"; do
	if ! [ -f "${i}.cntTable.gz" ]; then
		echo $i
		BAM_FILE="${BAMDIR}/${i}.bam"
		# nohup TEcount --GTF $gtf_genes --TE $gtf_transposable_elements --sortByPos --format BAM --mode multi --stranded reverse -b $BAM_FILE --project $i --outdir $OUTDIR &
		sbatch -n 4 --mem 15G -o slurm_${i}.out -e slurm_${i}.err --wrap "TEcount --GTF $gtf_genes --TE $gtf_transposable_elements --sortByPos --format BAM --mode multi --stranded no -b $BAM_FILE --project $i --outdir $OUTDIR"
	fi
done
