here::i_am("rna/TF2gene_coexpression/coexpression_TF_vs_gene_metacells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--TFs_file',             type="character",                               help='txt file with a list of TFs')
p$add_argument('--cor_test',         type="character",        help='Pearson or spearman')
p$add_argument('--remove_ExE_cells', action="store_true",                                 help='Remove ExE cells?')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args$sce <- file.path(io$basedir,"results/rna/metacells/all_cells/SingleCellExperiment_metacells.rds")
# args$metadata <- file.path(io$basedir,"results/rna/metacells/all_cells/metacells_metadata.txt.gz")
# args$TFs_file <- io$TFs_file
# args$cor_test <- "pearson"
# args$outdir <- file.path(io$basedir,"results/rna/coexpression")
# args$remove_ExE_cells <- TRUE
## END TEST ##

print(args)

# Sanity checks
stopifnot(args$cor_test%in%c("pearson","spearman"))

# I/O
dir.create(args$outdir, showWarnings=F)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

table(sample_metadata$celltype)

#########################
## Load RNA expression ##
#########################

# Load SingleCellExperiment
rna.sce <- readRDS(args$sce)[,sample_metadata$metacell]

##################################################
## Split RNA expression matrix into TF vs genes ##
##################################################

# TFs <- fread(args$TFs_file)[["gene"]] %>% unique
TFs <- fread(args$TFs_file)[[1]]# %>% unique

print("The following TFs are not found in the SingleCellExperiment:")
print(TFs[!TFs%in%toupper(rownames(rna.sce))])

rna_tfs.sce <- rna.sce[toupper(rownames(rna.sce))%in%TFs,]
rna_genes.sce <- rna.sce
rownames(rna_tfs.sce) <- toupper(rownames(rna_tfs.sce))

print(sprintf("Number of TFs: %s",nrow(rna_tfs.sce)))
print(sprintf("Number of genes: %s",nrow(rna_genes.sce)))

##########################
## Correlation analysis ##
##########################

if (args$cor_test=="pearson") {
  tf2gene_cor.mtx <- cor(t(as.matrix(logcounts(rna_tfs.sce))),t(as.matrix(logcounts(rna_genes.sce)))) %>% round(2)
  tf2tf_cor.mtx <- cor(t(as.matrix(logcounts(rna_tfs.sce))),t(as.matrix(logcounts(rna_tfs.sce)))) %>% round(2)
} else if (args$cor_test=="spearman") {
  if (is(logcounts(rna_tfs.sce),"sparseMatrix")) {
    tf2gene_cor.mtx <- SparseSpearmanCor2(t(logcounts(rna_tfs.sce)), t(logcounts(rna_genes.sce)))
    tf2tf_cor.mtx <- SparseSpearmanCor2(t(logcounts(rna_tfs.sce)), t(logcounts(rna_tfs.sce)))
  } else {
    tf2gene_cor.mtx <- cor(t(as.matrix(logcounts(rna_tfs.sce))), t(as.matrix(logcounts(rna_genes.sce))), method = "spearman")
    tf2tf_cor.mtx <- cor(t(as.matrix(logcounts(rna_tfs.sce))), t(as.matrix(logcounts(rna_tfs.sce))), method = "spearman")
  }
}

##########
## Save ##
##########

outfile <- sprintf("correlation_matrix_%s_tf2gene_metacells",args$cor_test)
outfile <- ifelse(args$remove_ExE_cells, paste0(outfile,"_no_ExE"),outfile)
outfile <- paste0(outfile,".rds")
saveRDS(tf2gene_cor.mtx, file.path(args$outdir,outfile))

outfile <- sprintf("correlation_matrix_%s_tf2tf_metacells",args$cor_test)
outfile <- ifelse(args$remove_ExE_cells, paste0(outfile,"_no_ExE"),outfile)
outfile <- paste0(outfile,".rds")
saveRDS(tf2tf_cor.mtx, file.path(args$outdir,outfile))