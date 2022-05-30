# Quantification of TF activities at the pseudobulk level
# https://www.archrproject.com/bookdown/chromvar-deviatons-enrichment-with-archr.html
# Note: this requires the previous execution of (...)/add_motif_annotation/archR_add_motif_annotation.R
here::i_am("atac/archR/chromvar/pseudobulk/run_chromvar_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(chromVAR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--peak_matrix_file',             type="character",           help='Precomputed peak matrix (pseudobulked)')
p$add_argument('--peak_metadata',             type="character",        help='Peak metadata file')
p$add_argument('--motif_annotation',             type="character",   help='Motif annotation')
p$add_argument('--motifmatcher',             type="character",     help='Motif annotation')
p$add_argument('--background_peaks',             type="character",        help='Background peaks')
p$add_argument('--min_peaks',             type="integer",  default=25,           help='Minimum number of peaks per TF')
p$add_argument('--outdir',          type="character",                    help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args$peak_matrix_file <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype_genotype/PeakMatrix/PeakMatrix_pseudobulk_with_replicates.rds")
# args$peak_matrix_file <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype_genotype/PeakMatrix/PeakMatrix_pseudobulk_with_replicates.rds")
# args$peak_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
# args$motif_annotation <- "CISBP"
# args$motifmatcher <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s-Scores.rds",args$motif_annotation))
# args$background_peaks <- file.path(io$basedir,"processed/atac/archR/Background-Peaks.rds")
# args$min_peaks <- 25
# args$outdir <- file.path(io$basedir,"results/atac/archR/chromvar/pseudobulk/test")
## END TEST ##

################################
## Load pseudobulk PeakMatrix ##
################################

print("Loading pseudobulk peak matrix...")

atac_peakMatrix_pseudobulk.se <- readRDS(args$peak_matrix_file)

# Load peak metadata
peak_metadata.dt <- fread(args$peak_metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

# Define peak names
# peak_names <- rowData(atac_peakMatrix_pseudobulk.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
# rownames(atac_peakMatrix_pseudobulk.se) <- peak_names

# Sanity checks
stopifnot(rownames(atac_peakMatrix_pseudobulk.se)%in%peak_metadata.dt$idx)

# Subset peaks
# atac_peakMatrix_pseudobulk.se <- atac_peakMatrix_pseudobulk.se[rownames(atac_peakMatrix_pseudobulk.se) %in% peak_metadata.dt$idx]
# peak_metadata.dt <- peak_metadata.dt %>%
#   .[idx%in%rownames(atac_peakMatrix_pseudobulk.se)] %>% setkey(idx) %>%
#   .[rownames(atac_peakMatrix_pseudobulk.se)]

# temporary
if (any(!c("start","strand")%in%colnames(rowData(atac_peakMatrix_pseudobulk.se)))) {
  tmp <- rownames(atac_peakMatrix_pseudobulk.se) %>% strsplit(":") %>% map_chr(2)
  rowData(atac_peakMatrix_pseudobulk.se)$start <- tmp %>% strsplit("-") %>% map_chr(1)
  rowData(atac_peakMatrix_pseudobulk.se)$end <- tmp %>% strsplit("-") %>% map_chr(2)
}

###############################
## Load motifmatcher results ##
###############################

print("Loading motifmatcher results...")

motifmatcher.se <- readRDS(args$motifmatcher)

###########################
## Load background peaks ##
###########################

print("Loading background peaks...")

bgdPeaks.se <- readRDS(args$background_peaks)
tmp <- rowRanges(bgdPeaks.se)
rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))

# stopifnot(rownames(atac_peakMatrix_pseudobulk.se)==rownames(bgdPeaks.se))
stopifnot(sort(rownames(atac_peakMatrix_pseudobulk.se))==sort(rownames(bgdPeaks.se)))

bgdPeaks.se <- bgdPeaks.se[rownames(atac_peakMatrix_pseudobulk.se),]

#########################################
## Use chromVAR default implementation ##
#########################################

print("Running chromVAR's default implementation...")

# prepare data for chromvar
assayNames(atac_peakMatrix_pseudobulk.se) <- "counts"

# Compute deviations
chromvar_deviations_chromvar.se <- computeDeviations(
  object = atac_peakMatrix_pseudobulk.se,
  annotations = motifmatcher.se,
  background_peaks = assay(bgdPeaks.se)
)

# Save
saveRDS(chromvar_deviations_chromvar.se, file.path(args$outdir,sprintf("chromVAR_deviations_%s_pseudobulk_chromvar.rds",args$motif_annotation)))

#########################################
## Use ArchR's chromVAR implementation ##
#########################################

print("Running ArchR's chromVAR implementation...")

source(here::here("atac/archR/chromvar/utils.R"))

featureDF <- data.frame(
  rowSums = rowSums(assay(atac_peakMatrix_pseudobulk.se)),
  start = rowData(atac_peakMatrix_pseudobulk.se)$start,
  end = rowData(atac_peakMatrix_pseudobulk.se)$end,
  GC = peak_metadata.dt$GC
)

chromvar_deviations_archr.se <- .customDeviations(
  countsMatrix = assay(atac_peakMatrix_pseudobulk.se),
  annotationsMatrix = as(assay(motifmatcher.se),"dgCMatrix"),
  backgroudPeaks = assay(bgdPeaks.se),
  expectation = featureDF$rowSums/sum(featureDF$rowSums),
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE
)

saveRDS(chromvar_deviations_archr.se, file.path(args$outdir,sprintf("chromVAR_deviations_%s_pseudobulk_archr.rds",args$motif_annotation)))

###################################################################
## Compare default chromVAR with ArchR's chromVAR implementation ##
###################################################################

# chromvar.deviations.dt_A <- assay(chromvar_deviations_chromvar.se,"z") %>% as.matrix %>% as.data.frame %>%
#   as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
#   melt(id.vars="motif", variable.name="celltype") %>%
#   .[,class:=factor("chromvar")]

# chromvar.deviations.dt_B <- assay(chromvar_deviations_archr.se) %>% as.matrix %>% as.data.frame %>%
#   as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
#   melt(id.vars="motif", variable.name="celltype") %>%
#   .[,class:=factor("archr")]

# to.plot <- rbind(chromvar.deviations.dt_A, chromvar.deviations.dt_B)

# opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype)]

# # motifs.to.plot <- unique(to.plot$motif)
# motifs.to.plot <- c("FOXA2","GATA2","TAL1")

# for (i in motifs.to.plot) {
  
#   to.plot.barplot <- to.plot[motif==i] %>% 
#     .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))]
  
#   p1 <- ggbarplot(to.plot.barplot, x="celltype", y="value", fill="celltype") +
#     scale_fill_manual(values=opts$celltype.colors) +
#     facet_wrap(~class, scales="free_y") +
#     geom_hline(yintercept=0, linetype="dashed") +
#     guides(x = guide_axis(angle = 90)) +
#     labs(x="", y=sprintf("%s chromVAR z-score",i)) +
#     theme(
#       legend.position = "none",
#       # axis.text.x = element_text(color="black", angle=40, hjust=1, size=rel(0.75)),
#       axis.text = element_text(color="black", size=rel(0.6)),
#       axis.title = element_text(color="black", size=rel(0.8))
#     )
  
#   to.plot.scatter <- to.plot.barplot %>% dcast(motif+celltype~class, value.var="value")
#   p2 <- ggscatter(to.plot.scatter, x="chromvar", y="archr", fill="celltype", shape=21, stroke=0.1, size=3.5, 
#                   add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#     stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
#     scale_fill_manual(values=opts$celltype.colors) +
#     labs(x="z-score (chromVAR)", y="z-score (ArchR)") +
#     theme(
#       legend.position = "none",
#       axis.text = element_text(size=rel(0.7)),
#       axis.title = element_text(size=rel(0.85))
#     )
  
#   p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1, rel_widths = c(2/3,1/3))
  
#   pdf(file.path(args$outdir,sprintf("%s_chromvar_singlecell_vs_pseudobulk.pdf",i)), width=13, height=5)
#   print(p)
#   dev.off()
# }

################
## START TEST ##
################

# library(BiocParallel)
# R.utils::sourceDirectory("/Users/ricard/git/chromVAR/R/", verbose=T, modifiedOnly=FALSE)
# R.utils::sourceDirectory("/homes/ricard/git/chromVAR/R/", verbose=T, modifiedOnly=FALSE)

# object <- atac_peakMatrix_pseudobulk.se
# background_peaks = assay(bgdPeaks.se)
# annotations = motifmatcher.se
# rowData = NULL
# colData = NULL
# 
# counts_mat <- counts(atac_peakMatrix_pseudobulk.se)
# expectation = computeExpectations(counts_mat)
# peak_indices <- convert_to_ix_list(annotationMatches(annotations))
# 
# peak_set = peak_indices[[2]]
