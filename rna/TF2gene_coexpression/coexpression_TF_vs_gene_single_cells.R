here::i_am("rna/TF2gene_coexpression/coexpression_TF_vs_gene_single_cells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',              type="character",        help='SingleCellExperiment file')
p$add_argument('--metadata',         type="character",        help='Cell metadata file')
p$add_argument('--TFs_file',         type="character",        help='txt file with a list of TFs')
p$add_argument('--remove_ExE_cells', action="store_true",     help='Remove ExE cells?')
p$add_argument('--denoise',          action="store_true",     help='Apply KNN denoise?')
p$add_argument('--knn',              type="integer",          help='Number of KNN for denoise')
p$add_argument('--cor_test',         type="character",        help='Pearson or spearman')
p$add_argument('--outdir',          type="character",         help='Output file')
p$add_argument('--test_mode',       action="store_true",      help='Test mode? subset data')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$TFs_file <- "/Users/argelagr/data/mm10_regulation/TFs/TFs.txt" # "/bi/group/reik/ricard/data/mm10_regulation/TFs/TFs.txt"
# args$outdir <- file.path(io$basedir,"results/rna/coexpression")
# args$denoise <- FALSE
# args$remove_ExE_cells <- TRUE
# args$cor_test <- "pearson"
# args$knn <- 50
# args$test_mode <- TRUE
## END TEST ##

print(args)

# Sanity checks
stopifnot(args$cor_test%in%c("pearson","spearman"))

# I/O
dir.create(args$outdir, showWarnings=F)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE]

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>% .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

if (args$test_mode) {
	print("Test mode activated...")
	sample_metadata <- sample_metadata %>% head(n=150)
}

table(sample_metadata$celltype)

#########################
## Load RNA expression ##
#########################

# Load SingleCellExperiment
rna.sce <- load_SingleCellExperiment(
	file = args$sce, 
	cells = sample_metadata$cell, 
	remove_non_expressed_genes = FALSE, 
	normalise = TRUE
)

# Update colData
tmp <- sample_metadata %>% .[cell%in%colnames(rna.sce)] %>% setkey(cell) %>% .[colnames(rna.sce)]
stopifnot(tmp$cell == colnames(rna.sce))
colData(rna.sce) <- tmp %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(rna.sce),] %>% DataFrame()

#############
## Denoise ##
#############

if (args$denoise) {
  
	print("denoise...")

	# Feature selection and dimensionality reduction  
	decomp <- modelGeneVar(rna.sce)
	hvgs <- decomp[order(decomp$FDR),] %>% head(n=2500) %>% rownames
	rna_filt.sce <- runPCA(rna.sce[hvgs,], ncomponents = 50, ntop=2500)

	# kNN denoise
	logcounts(rna.sce) <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(rna.sce)), D=pdist(reducedDim(rna_filt.sce,"PCA")), k=args$knn)
}

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

outfile <- sprintf("correlation_matrix_%s_tf2gene_single_cells",args$cor_test)
outfile <- ifelse(args$remove_ExE_cells, paste0(outfile,"_no_ExE"),outfile)
outfile <- ifelse(args$denoise, paste0(outfile,"_denoised"),outfile)
outfile <- paste0(outfile,".rds")
saveRDS(tf2gene_cor.mtx, file.path(args$outdir,outfile))

outfile <- sprintf("correlation_matrix_%s_tf2tf_single_cells",args$cor_test)
outfile <- ifelse(args$remove_ExE_cells, paste0(outfile,"_no_ExE"),outfile)
outfile <- ifelse(args$denoise, paste0(outfile,"_denoised"),outfile)
outfile <- paste0(outfile,".rds")
saveRDS(tf2tf_cor.mtx, file.path(args$outdir,outfile))

##########
## TEST ##
##########

# tf2gene_cor.mtx <- readRDS(file.path(args$outdir,"correlation_matrix_tf2gene.rds"))

# i <- "FOXA2"
# j <- "Cab39l"

# to.plot <- data.table(
#   TF = logcounts(rna_tfs.sce[i,])[1,],
#   target_gene = logcounts(rna_genes.sce[j,])[1,],
#   celltype = colnames(rna_tfs.sce)
# )


# ggscatter(to.plot, x="TF", y="target_gene", fill="celltype", size=4, shape=21, 
#           add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#   stat_cor(method = "pearson") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   labs(x=sprintf("%s expression",i), y=sprintf("%s expression",j)) +
#   guides(fill=F) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size=rel(0.85)),
#     axis.text = element_text(size=rel(0.7))
#   )


##########
## TEST ##
##########

# system.time(foo <- SparseSpearmanCor2(t(logcounts(rna_tfs.sce)), t(logcounts(rna_genes.sce))))
# system.time(bar <- cor(t(as.matrix(logcounts(rna_tfs.sce))), t(as.matrix(logcounts(rna_genes.sce))), method = "spearman"))
# 
# dimnames(foo) <- dimnames(bar)
# foo[1:3,1:3]
# bar[1:3,1:3]
# foo[is.na(foo)] <- 0
# bar[is.na(bar)] <- 0
# norm(foo - bar, type = "2")
