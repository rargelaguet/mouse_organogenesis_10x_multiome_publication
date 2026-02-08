##############
## Settings ##
##############

# BAM_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac"
# samples=( "pA_Control_REP1" "pA_Control_REP2" "pA_Control_REP3" "pA_TG_REP1" "pA_TG_REP2" "pA_TG_REP3" "pB_Control_REP1" "pB_Control_REP2" "pB_Control_REP3" "pB_TG_REP1" "pB_TG_REP2" "pB_TG_REP3" "pC_Control_REP1" "pC_Control_REP2" "pC_Control_REP3" "pC_TG_REP1" "pC_TG_REP2" "pC_TG_REP3" "pD_Control_REP1" "pD_Control_REP2" "pD_Control_REP3" "pD_TG_REP1" "pD_TG_REP2" "pD_TG_REP3" "pE_Control_REP1" "pE_Control_REP2" "pE_TG_REP1" "pE_TG_REP2" "pF_Control_REP1" "pF_Control_REP2" "pF_Control_REP3" "pF_TG_REP1" "pF_TG_REP2" "pF_TG_REP3" "pG_Control_REP1" "pG_Control_REP2" "pG_Control_REP3" "pG_TG_REP1" "pG_TG_REP2" )

BAM_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac/merged"
samples=( "Control_p11" "Control_p14" "Control_p19" "Control_p22" "Control_p24" "Control_p29" "Control_p8" "TG_p11" "TG_p14" "TG_p19" "TG_p22" "TG_p24" "TG_p29" "TG_p8" )

###########
## MACS2 ##
###########

# 1) bedtools bamtobed -i $bam > ${name}_pe.bed 
# 2) macs2 callpeak -t ${i}_pe.bed -n ${i}_narrow -f BED -g ${species} -q 0.01 --nomodel --shift -75 --extsize 150 --call-summits --keep-dup all

# i="pD_TG_REP1"
for i in "${samples[@]}"; do
    if ! [ -f "${OUTDIR}/${i}.macs2_summits.bed" ]; then
        echo $i
        # nohup (...) &> ${i}_nohup.out &
        sbatch -n 4 --mem 15G -o slurm_${i}.out --wrap "(...)"
    fi
done


###########
## MACS3 ##
###########

OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/results/atac/peak_calling/macs3"

# i="pD_TG_REP1"
for i in "${samples[@]}"; do
    if ! [ -f "${OUTDIR}/${i}.macs3_summits.bed" ]; then
        echo $i
        nohup macs3 callpeak -f BAMPE --call-summits -t ${BAM_DIR}/${i}.bam -n ${i} --outdir ${OUTDIR} -B -q 0.01 &> ${i}_nohup.out &
        # sbatch -n 4 --mem 15G -o slurm_${i}.out --wrap "macs3 callpeak -f BAMPE --call-summits -t ${BAM_DIR}/${i}.bam -n ${i}.macs3 --outdir ${OUTDIR} -B -q 0.01"
    fi
done

# filter by q-value and extract bed files
for i in "${samples[@]}"; do
    awk -v OFS="\t" '$5>=50 {print}' ${i}_peaks.narrowPeak | cut -f 1,2,3 > ${i}_peaks_filt.bed
done

# merge peaks
bedops -m nk_peaks.bed nk_stim_peaks.bed > nk_merged_peaks.bed

bedops -m Control_p11_peaks_filt.bed Control_p14_peaks_filt.bed Control_p19_peaks_filt.bed Control_p22_peaks_filt.bed Control_p24_peaks_filt.bed Control_p29_peaks_filt.bed Control_p8_peaks_filt.bed TG_p11_peaks_filt.bed TG_p14_peaks_filt.bed TG_p19_peaks_filt.bed TG_p22_peaks_filt.bed TG_p24_peaks_filt.bed TG_p29_peaks_filt.bed TG_p8_peaks_filt.bed > merged_peaks.bed

####################
## MACS3 HMMRATAC ##
####################

OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/results/atac/peak_calling/macs3_hmmratac"

# i="pD_TG_REP1"
for i in "${samples[@]}"; do
    if ! [ -f "${OUTDIR}/${i}.macs3_summits.bed" ]; then
        echo $i
        # nohup macs3 hmmratac -i ${BAM_DIR}/${i}.bam -n ${i}.macs3_hmmratac &> ${i}_nohup.out &
        sbatch -n 4 --mem 15G -o slurm_${i}.out --wrap "macs3 hmmratac -i ${BAM_DIR}/${i}.bam -n ${i}.macs3_hmmratac"
    fi
done

# Filter peaks
# awk -v OFS="\t" '$13>=10 {print}' NAME_peaks.gappedPeak > NAME.filteredPeaks.gappedPeak

# create bed file
