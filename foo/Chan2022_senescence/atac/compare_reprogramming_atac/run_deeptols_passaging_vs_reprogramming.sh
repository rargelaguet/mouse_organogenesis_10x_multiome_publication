conda activate genomic_ops # local
# conda activate ucsc # dockyard

#########################################################################
## Quantify ATAC-seq signal of BJ reprogramming using BJ passage peaks ##
#########################################################################

# BED_FILE="/share/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/archR/peaks_archR_macs2.bed"
BED_FILE="/Users/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/jack/peak_calls.bed"
BIGWIG_PREFIX="/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/archR/GroupBigWigs/cluster/symbolic"
OUTDIR="/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/deeptools_pseudobulk"

mkdir -p $OUTDIR; cd $OUTDIR

OUTPREFIX="${OUTDIR}/multiBigwigSummary_BJ_passaging"

nohup multiBigwigSummary BED-file \
       -b $BIGWIG_PREFIX/*.bw \
       --BED $BED_FILE \
       -o ${OUTPREFIX}.npz \
       --outRawCounts ${OUTPREFIX}.tab \
       -p 2 &

###############################################################################################
## Quantify ATAC-seq signal of BJ passaging using BJ reprogramming peaks (from BigWig files) ##
###############################################################################################

# BED_FILE="/share/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/archR/peaks_archR_macs2.bed"
BED_FILE="/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/archR/PeakCalls/peaks_archR_macs2.bed"
BIGWIG_PREFIX="/Users/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/bigwig"
OUTDIR="/Users/rargelaguet/data/BJ_passaging_CHAMP/results/atac/deeptools"

mkdir -p $OUTDIR; cd $OUTDIR

OUTPREFIX="${OUTDIR}/multiBigwigSummary_BJ_reprogramming"

nohup multiBigwigSummary BED-file \
       -b $BIGWIG_PREFIX/*.bw \
       --BED $BED_FILE \
       -o ${OUTPREFIX}.npz \
       --outRawCounts ${OUTPREFIX}.tab \
       -p 1 &


############################################################################################
## Quantify ATAC-seq signal of BJ passaging using BJ reprogramming peaks (from BAM files) ##
############################################################################################

# # BED_FILE="/share/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/archR/peaks_archR_macs2.bed"
# BED_FILE="/share/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/archR/peaks_archR_macs2_no_chr.bed"
# BAM_DIR="/share/rargelaguet/data/BJ_passaging_CHAMP/bam/atac/merged"
# OUTDIR="/share/rargelaguet/data/BJ_passaging_CHAMP/results/atac/compare_reprogramming"

# samples=( "Control_p11" "Control_p14" "Control_p19" "Control_p22" "Control_p24" "Control_p29" "Control_p8" "TG_p11" "TG_p14" "TG_p19" "TG_p22" "TG_p24" "TG_p29" "TG_p8" )

# mkdir -p $OUTDIR; cd $OUTDIR

# # Create an array of bam files
# bam_files=()
# for sample in "${samples[@]}"; do
#   bam_files+=("${BAM_DIR}/${sample}.bam")
# done

# nohup bedtools multicov -bams ${bam_files[@]} -bed $BED_FILE > ${OUTDIR}/atac_coverage_peaks.count &
