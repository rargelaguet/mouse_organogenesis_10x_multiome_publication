
# BAM="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac/pA_TG_REP3.bam"
BAM_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac/merged"
PEAKS="/share/rargelaguet/data/BJ_passaging_CHAMP/results/atac/peak_calling/macs3/merged_peaks.bed"
BLACKLIST="/share/rargelaguet/data/hg38_sequence/GRCh38_blacklist.bed"
GENOME="/share/rargelaguet/data/hg38_sequence/whole_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias"
CORES=8

mkdir -p $OUTDIR; cd $OUTDIR

# samples=( "p8_Control_REP1" "p8_Control_REP2" "p8_Control_REP3" "p8_TG_REP1" "p8_TG_REP2" "p8_TG_REP3" "p11_Control_REP1" "p11_Control_REP2" "p11_Control_REP3" "p11_TG_REP1" "p11_TG_REP2" "p11_TG_REP3" "p14_Control_REP1" "p14_Control_REP2" "p14_Control_REP3" "p14_TG_REP1" "p14_TG_REP2" "p14_TG_REP3" "p19_Control_REP1" "p19_Control_REP2" "p19_Control_REP3" "p19_TG_REP1" "p19_TG_REP2" "p19_TG_REP3" "p22_Control_REP1" "p22_Control_REP2" "p22_Control_REP3" "p22_TG_REP1" "p22_TG_REP2" "p22_TG_REP3" "p24_Control_REP1" "p24_Control_REP2" "p24_Control_REP3" "p24_TG_REP1" "p24_TG_REP2" "p24_TG_REP3" "p29_Control_REP1" "p29_Control_REP2" "p29_Control_REP3" "p29_TG_REP1" "p29_TG_REP2" "p29_TG_REP3" )
# samples=( "p8_Control" "p11_Control" "p14_Control" "p19_Control" "p22_Control" "p24_Control" "p29_Control" "p8_TG" "p11_TG" "p14_TG" "p19_TG" "p22_TG" "p24_TG" "p29_TG" )
samples=( "p8_Control" "p8_TG" )

#######################
## TOBIAS ATACorrect ##
#######################

# i="p8_Control"
for i in "${samples[@]}"; do
    if ! [ -f "${OUTDIR}/${i}_corrected.bw" ]; then
        echo $i
        nohup TOBIAS ATACorrect --bam ${BAM_DIR}/${i}.bam --genome $GENOME --peaks $PEAKS --blacklist $BLACKLIST --outdir $OUTDIR --cores $CORES &> ${i}_nohup.out &
        # sbatch -n 4 --mem 30G -o slurm_${i}.out --wrap "Rscript $SCRIPT --sample $i"
    fi
done


############################
## TOBIAS FootprintScores ##
############################

# i="p8_Control"
for i in "${samples[@]}"; do
    if ! [ -f "${OUTDIR}/${i}_footprints.bw" ]; then
    	BW_CORRECTED=
        echo $i
        nohup TOBIAS FootprintScores --signal ${OUTDIR}/${i}_corrected.bw --regions $PEAKS --output ${OUTDIR}/${i}_footprints.bw --cores 8 &> ${i}_FootprintScores_nohup.out &
        # sbatch -n 4 --mem 30G -o slurm_${i}.out --wrap "Rscript $SCRIPT --sample $i"
    fi
done


######################
## TOBIAS BINDetect ##
######################

FOOTPRINTS="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias/pA_TG_REP3_footprints.bw"
MOTIFS="/share/rargelaguet/data/hg38_regulation/tf_motifs/jaspar/vertebrates/all_motifs.jaspar"
OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias/BINDetect"

mkdir $OUTDIR; cd $OUTDIR

# TOBIAS BINDetect --motifs $MOTIFS --signals $FOOTPRINTS --genome $GENOME --peaks $PEAKS --outdir $OUTDIR --cond_names Bcell Tcell --cores 8
nohup TOBIAS BINDetect --motifs $MOTIFS --signals $FOOTPRINTS --genome $GENOME --peaks $PEAKS --outdir $OUTDIR --cores $CORES &


#####################################
## TOBIAS BINDetect, Control vs Tg ##
#####################################

samples=( "p8_Control" "p8_TG" )

BASEDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias"
MOTIFS="/share/rargelaguet/data/hg38_regulation/tf_motifs/jaspar/vertebrates/all_motifs.jaspar"
OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias/BINDetect"

mkdir $OUTDIR; cd $OUTDIR

nohup TOBIAS BINDetect --motifs $MOTIFS --signals $BASEDIR/p8_Control_footprints.bw $BASEDIR/p8_TG_footprints.bw --genome $GENOME --peaks $PEAKS --outdir $OUTDIR --cond_names Control TG --cores 8 &> Control-vs_TG_BINDetect_nohup.out &


###################
## PlotAggregate ##
###################

BASEDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias"
BOUND_BED="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias/BINDetect/ATF4_MA0833.2/bed/ATF4_MA0833.2_TG_bound.bed"
UNBOUND_BED="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias/BINDetect/ATF4_MA0833.2/bed/ATF4_MA0833.2_TG_unbound.bed"
FOOTPRINTS="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias/pA_TG_REP3_footprints.bw"
OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/footprinting_tobias/BINDetect"


TOBIAS PlotAggregate --TFBS test_data/BATF_all.bed  --signals test_data/Bcell_corrected.bw test_data/Tcell_corrected.bw --output BATFJUN_footprint_comparison_all.pdf --share_y both --plot_boundaries --signal-on-x
TOBIAS PlotAggregate --TFBS $BOUND_BED $UNBOUND_BED --signals test_data/Bcell_uncorrected.bw test_data/Bcell_expected.bw test_data/Bcell_corrected.bw --output IRF1_footprint.png --share_y sites --plot_boundaries
