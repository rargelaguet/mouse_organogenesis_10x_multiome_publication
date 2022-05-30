here::i_am("rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/run_cor_TFexpr_vs_peakAcc_metacells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (metacells)') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (metacells)') 
p$add_argument('--motif2gene',  type="character",              help='Motif annotation') 
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation') 
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--test',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/metacells/SingleCellExperiment_metacells.rds")
# args$atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/metacells/PeakMatrix_summarized_experiment_metacells.rds")
# args$sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")
# args$motif_annotation <- "CISBP"
# args$motifmatcher <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,args$motif_annotation)
# args$motif2gene <- sprintf("%s/Annotations/%s_TFs.txt.gz",io$archR.directory,args$motif_annotation)
# args$test <- FALSE
# args$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/trajectories/nmp/test")
## END TEST ##

## START TEST ##
# args <- list()
# io$basedir <- file.path(io$basedir,"test")
# args$sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")
# args$atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/metacells/trajectories/nmp/PeakMatrix/PeakMatrix_summarized_experiment_metacells.rds")
# args$sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")
# args$motif_annotation <- "CISBP"
# args$motifmatcher <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s-Scores.rds",args$motif_annotation))
# args$motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s_motif2gene.txt.gz",args$motif_annotation))
# args$test <- FALSE
# args$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/trajectories/nmp/test")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F, recursive = T)

#######################
## Load RNA and ATAC ##
#######################

# Load SingleCellExperiment
rna_metacells.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_peakMatrix_metacells.se <- readRDS(args$atac_peak_matrix)

# Normalise ATAC data
assayNames(atac_peakMatrix_metacells.se) <- "counts"
assay(atac_peakMatrix_metacells.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_metacells.se),2,colSums(assay(atac_peakMatrix_metacells.se),na.rm=T),"/"))+1)

# Make sure that samples are consistent
samples <- intersect(colnames(rna_metacells.sce),colnames(atac_peakMatrix_metacells.se))
rna_metacells.sce <- rna_metacells.sce[,samples]
atac_peakMatrix_metacells.se <- atac_peakMatrix_metacells.se[,samples]

print(sprintf("Number of metacells: %s",length(samples)))

###############################
## Load motifmatcher results ##
###############################

motifmatcher.se <- readRDS(args$motifmatcher)

# Subset peaks
stopifnot(sort(rownames(motifmatcher.se))==sort(rownames(atac_peakMatrix_metacells.se)))
motifmatcher.se <- motifmatcher.se[rownames(atac_peakMatrix_metacells.se),]

###########################
## Load motif annotation ##
###########################

motif2gene.dt <- fread(args$motif2gene)

################
## Filter TFs ##
################

motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
motifmatcher.se <- motifmatcher.se[,motifs]
motif2gene.dt <- motif2gene.dt[motif%in%motifs]

genes <- intersect(toupper(rownames(rna_metacells.sce)),motif2gene.dt$gene)
rna_tf_metacells.sce <- rna_metacells.sce[str_to_title(genes),]
rownames(rna_tf_metacells.sce) <- toupper(rownames(rna_tf_metacells.sce))
motif2gene.dt <- motif2gene.dt[gene%in%genes]

# Manually remove some motifs
motif2gene.dt <- motif2gene.dt[!motif%in%c("T_789")]

# Remove duplicated gene-motif pairs
genes.to.remove <- names(which(table(motif2gene.dt$gene)>1))
cat(sprintf("Removing %d TFs that have duplicated gene-motif pairs:\n%s", length(genes.to.remove), paste(genes.to.remove, collapse=", ")))
motif2gene.dt <- motif2gene.dt[!gene%in%genes.to.remove]
rna_tf_metacells.sce <- rna_tf_metacells.sce[rownames(rna_tf_metacells.sce)%in%motif2gene.dt$gene]
motifmatcher.se <- motifmatcher.se[,colnames(motifmatcher.se)%in%motif2gene.dt$motif]
stopifnot(table(motif2gene.dt$gene)==1)

#########################################################
## Correlate peak accessibility with TF RNA expression ##
#########################################################

# Sanity checks
stopifnot(colnames(rna_tf_metacells.sce)==colnames(atac_peakMatrix_metacells.se))

TFs <- rownames(rna_tf_metacells.sce)
if (args$test) {
  TFs <- TFs %>% head(n=5)
}

# Prepare output data objects
cor.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_peakMatrix_metacells.se), ncol=length(TFs))
pvalue.mtx <- matrix(as.numeric(NA), nrow=nrow(atac_peakMatrix_metacells.se), ncol=length(TFs))
rownames(cor.mtx) <- rownames(atac_peakMatrix_metacells.se); colnames(cor.mtx) <- TFs
dimnames(pvalue.mtx) <- dimnames(cor.mtx)

# i <- "GATA1"
for (i in TFs) {
  motif_i <- motif2gene.dt[gene==i,motif]
  
  print(i)
  all_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,motif_i],"motifMatches")==1)]
  
  # calculate correlations
  corr_output <- psych::corr.test(
    x = t(logcounts(rna_tf_metacells.sce[i,])), 
    y = t(assay(atac_peakMatrix_metacells.se[all_peaks_i,],"logcounts")), 
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
  rowData = rowData(atac_peakMatrix_metacells.se)
)
saveRDS(to.save, file.path(args$outdir,sprintf("%s_cor_TFexpr_vs_peakAcc.rds",args$motif_annotation)))

##########
## TEST ##
##########

# assay(motifmatcher.se["chr5:4894870-4895470","TBXT"],"motifMatches")
# corr_output$r[1,]["chr5:4894870-4895470"]
# cor.mtx["chr5:4894870-4895470","T"]
