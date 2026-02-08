conda activate genomic_ops # local
# conda activate ucsc # dockyard

#####################
## Define settings ##
#####################

BED_FILE="/Users/rargelaguet/data/hg38_regulation/chromHMM/HDF_6_dense.merged.bed"
BIGWIG_PREFIX="/Users/rargelaguet/data/BJ_passaging_CHAMP/processed/atac/bigwig"
OUTDIR="/Users/rargelaguet/data/BJ_passaging_CHAMP/results/atac/deeptools"

# BED_FILE="/share/rargelaguet/data/hg38_regulation/chromHMM/HDF_6_dense.merged.bed"
# # BED_FILE="/share/rargelaguet/data/hg38_regulation/transposable_elements/ucsc_repeatmasker/bed/name/LTR7Y.bed"
# BIGWIG_PREFIX="/share/rargelaguet/data/single_cell_mptr/processed/atac/archR/GroupBigWigs/sample_day_celltype"
# OUTDIR="/share/rargelaguet/data/single_cell_mptr/results/atac/archR/deeptools_pseudobulk/unnormalised"

THREADS=1

mkdir -p $OUTDIR; cd $OUTDIR

###################
## computeMatrix ##
###################

OUTPREFIX="${OUTDIR}/computeMatrix_chromHMM"

nohup computeMatrix reference-point \
       -b 500 -a 500 \
       -R $BED_FILE \
       -S $BIGWIG_PREFIX/*.bw \
       -o ${OUTPREFIX}.gz \
       --outFileSortedRegions ${OUTPREFIX}.tab \
       -p $THREADS &


#################################
## multiBigwigSummary BED-file ##
#################################

OUTPREFIX="${OUTDIR}/multiBigwigSummary_chromHMM"

nohup multiBigwigSummary BED-file \
       -b $BIGWIG_PREFIX/*.bw \
       --BED $BED_FILE \
       -o ${OUTPREFIX}.npz \
       --outRawCounts ${OUTPREFIX}.tab \
       -p $THREADS &

############################
## multiBigwigSummary bin ##
############################

# multiBigwigSummary bins --binSize 100 -p 1 -b $BIGWIG_PREFIX/*.bw -o test_multiBigwigSummary.npz 
# multiBigwigSummary bins --binSize 100 -p 1 -b 25F_gapdh_high.bw 25F_gapdh_low.bw 28F_gapdh_high.bw 28F_gapdh_low.bw 52F_gapdh_high.bw -o results.npz 