here::i_am("rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/run_cor_TFexpr_vs_peakAcc_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--motif2gene',  type="character",              help='Motif2gene') 
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation') 
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--test',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
# args$atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype/PeakMatrix/pseudobulk_PeakMatrix_summarized_experiment.rds")
# args$motif_annotation <- "CISBP"
# args$motifmatcher <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s-Scores.rds",args$motif_annotation))
# args$motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s_motif2gene.txt.gz",args$motif_annotation))
# args$test <- FALSE
# args$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/test")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F, recursive = T)

#######################
## Load RNA and ATAC ##
#######################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_peakMatrix_pseudobulk.se <- readRDS(args$atac_peak_matrix)

# Normalise ATAC data
assayNames(atac_peakMatrix_pseudobulk.se) <- "counts"
assay(atac_peakMatrix_pseudobulk.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_pseudobulk.se),2,colSums(assay(atac_peakMatrix_pseudobulk.se),na.rm=T),"/"))+1)

# Make sure that samples are consistent
samples <- intersect(colnames(rna_pseudobulk.sce),colnames(atac_peakMatrix_pseudobulk.se))
rna_pseudobulk.sce <- rna_pseudobulk.sce[,samples]
atac_peakMatrix_pseudobulk.se <- atac_peakMatrix_pseudobulk.se[,samples]

###############################
## Load motifmatcher results ##
###############################

motifmatcher.se <- readRDS(args$motifmatcher)

# Subset peaks
stopifnot(sort(rownames(motifmatcher.se))==sort(rownames(atac_peakMatrix_pseudobulk.se)))
motifmatcher.se <- motifmatcher.se[rownames(atac_peakMatrix_pseudobulk.se),]

################################
## Load motif2gene annotation ##
################################

motif2gene.dt <- fread(args$motif2gene)

################
## Filter TFs ##
################

motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
motifmatcher.se <- motifmatcher.se[,motifs]
motif2gene.dt <- motif2gene.dt[motif%in%motifs]

genes <- intersect(toupper(rownames(rna_pseudobulk.sce)),motif2gene.dt$gene)
rna_tf_pseudobulk.sce <- rna_pseudobulk.sce[str_to_title(genes),]
rownames(rna_tf_pseudobulk.sce) <- toupper(rownames(rna_tf_pseudobulk.sce))
motif2gene.dt <- motif2gene.dt[gene%in%genes]

# Manually remove some motifs
motif2gene.dt <- motif2gene.dt[!motif%in%c("T_789")]

# Remove duplicated gene-motif pairs
genes.to.remove <- names(which(table(motif2gene.dt$gene)>1))
cat(sprintf("Removing %d TFs that have duplicated gene-motif pairs:\n%s", length(genes.to.remove), paste(genes.to.remove, collapse=", ")))
motif2gene.dt <- motif2gene.dt[!gene%in%genes.to.remove]
rna_tf_pseudobulk.sce <- rna_tf_pseudobulk.sce[rownames(rna_tf_pseudobulk.sce)%in%motif2gene.dt$gene]
motifmatcher.se <- motifmatcher.se[,colnames(motifmatcher.se)%in%motif2gene.dt$motif]
stopifnot(table(motif2gene.dt$gene)==1)

#########################################################
## Correlate peak accessibility with TF RNA expression ##
#########################################################

# Sanity checks
stopifnot(colnames(rna_tf_pseudobulk.sce)==colnames(atac_peakMatrix_pseudobulk.se))

TFs <- rownames(rna_tf_pseudobulk.sce)
if (args$test) {
  TFs <- TFs %>% head(n=5)
}

# Prepare output data objects
cor.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_peakMatrix_pseudobulk.se), ncol=length(TFs))
pvalue.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_peakMatrix_pseudobulk.se), ncol=length(TFs))
rownames(cor.mtx) <- rownames(atac_peakMatrix_pseudobulk.se); colnames(cor.mtx) <- TFs
dimnames(pvalue.mtx) <- dimnames(cor.mtx)

# i <- "GATA1"
for (i in TFs) {
  motif_i <- motif2gene.dt[gene==i,motif]
  
  print(i)
  all_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,motif_i],"motifMatches")==1)]
  
  # calculate correlations
  corr_output <- psych::corr.test(
    x = t(logcounts(rna_tf_pseudobulk.sce[i,])), 
    y = t(assay(atac_peakMatrix_pseudobulk.se[all_peaks_i,],"logcounts")), 
    ci = FALSE
  )
  
  # Fill matrices
  cor.mtx[all_peaks_i,i] <- round(corr_output$r[1,],3)
  pvalue.mtx[all_peaks_i,i] <- round(corr_output$p[1,],10)
}

##########
## Save ##
##########

to.save <- SummarizedExperiment(
  assays = SimpleList("cor" = dropNA(cor.mtx), "pvalue" = dropNA(pvalue.mtx)),
  rowData = rowData(atac_peakMatrix_pseudobulk.se)
)
saveRDS(to.save, file.path(args$outdir,sprintf("%s_cor_TFexpr_vs_peakAcc.rds",args$motif_annotation)))

##########
## TEST ##
##########

# assay(motifmatcher.se["chr15:36706260-36706860","TBXT_121"],"motifScores")
# assay(motifmatcher.se["chr16:11469994-11470594","TBXT_121"],"motifCounts")

