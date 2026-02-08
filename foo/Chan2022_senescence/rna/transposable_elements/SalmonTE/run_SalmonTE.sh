
##############
## Settings ##
##############

SALMON="/share/rargelaguet/git/SalmonTE/salmon/linux/bin/salmon"
FASTQ_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/fastq/symbolic"
OUTPUT_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/SalmonTE"
GENOME_INDEX="/share/rargelaguet/git/SalmonTE/reference/hs"

mkdir $OUTPUT_DIR; cd $OUTPUT_DIR

samples=( "p8_Control_REP1" "p8_Control_REP2" "p8_Control_REP3" "p8_TG_REP1" "p8_TG_REP2" "p8_TG_REP3" "p11_Control_REP1" "p11_Control_REP2" "p11_Control_REP3" "p11_TG_REP1" "p11_TG_REP2" "p11_TG_REP3" "p14_Control_REP1" "p14_Control_REP2" "p14_Control_REP3" "p14_TG_REP1" "p14_TG_REP2" "p14_TG_REP3" "p19_Control_REP1" "p19_Control_REP2" "p19_Control_REP3" "p19_TG_REP1" "p19_TG_REP2" "p19_TG_REP3" "p22_Control_REP1" "p22_Control_REP2" "p22_Control_REP3" "p22_TG_REP1" "p22_TG_REP2" "p22_TG_REP3" "p24_Control_REP1" "p24_Control_REP2" "p24_Control_REP3" "p24_TG_REP1" "p24_TG_REP2" "p24_TG_REP3" "p29_Control_REP1" "p29_Control_REP2" "p29_Control_REP3" "p29_TG_REP1" "p29_TG_REP2" "p29_TG_REP3" )

##################
## Run SalmonTE ##
##################

# i="p8_Control_REP1"
for i in "${samples[@]}"; do
	if ! [ -f "${OUTPUT_DIR}/${i}/quant.sf" ]; then
		# echo $i
    	# nohup $SALMON quant -q -i $GENOME_INDEX -l A -r ${FASTQ_DIR}/${i}_R1.fastq.gz -o ${OUTPUT_DIR}/${i} -p 1 &
    	sbatch -n 4 --mem 15G -o slurm_${i}.out -e slurm_${i}.err --wrap "$SALMON quant -q -i $GENOME_INDEX -l A -r ${FASTQ_DIR}/${i}_R1.fastq.gz -o ${OUTPUT_DIR}/${i} -p 1"
	fi
done

