# Cell type assignments are done by mapping the RNA expression profiles to the PijuanSala2019 atlas.
# there are some cells which did not pass QC for RNA expression but did pass QC on the ATAC modality
# This script predicts the cell type label for these cells by using a kNN prediction

# Note that this should be applied only to high-quality ATAC cells.

here::i_am("atac/archR/celltype_assignment/archR_celltype_assignment.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

library(irlba)

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--atac_peak_matrix',             type="character",            help='Atac peak matrix')
p$add_argument('--input_celltype_column',             type="character",             help='Input column for the celltype information')
p$add_argument('--output_celltype_column',             type="character",            help='Output column for the celltype information')
p$add_argument('--atac_feature_stats',             type="character",            help='')
p$add_argument('--k',             type="integer",  default=25,                      help='Number of kNN for the celltype prediction')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- file.path(io$basedir,"test/results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$atac_peak_matrix <- file.path(io$basedir,"test/processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
# args$input_celltype_column <- "celltype.mapped"
# args$output_celltype_column <- "celltype.predicted"
# args$atac_feature_stats <- file.path(io$basedir, "test/results/atac/archR/feature_stats/PeakMatrix_celltype.mapped_stats.txt.gz")
# args$k <- 25
# args$outdir <- file.path(io$basedir,"test/results/atac/archR/celltype_assignment")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

# Options
opts$nfeatures <- 50000
opts$ndims <- 100
opts$scale_dims <- TRUE
opts$remove_dim_cor_seq_depth <- TRUE

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata)

stopifnot(args$input_celltype_column%in%colnames(sample_metadata))

# Select high-quality ATAC cells that passed QC for RNA expression
sample_metadata <- sample_metadata %>%
  .[pass_atacQC==TRUE] %>% 
  # .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE & nFrags_atac>=1e4] %>% 
  setnames(args$input_celltype_column,"celltype")

if (sum(is.na(sample_metadata$celltype))==0) {
  stop("No cells are missing cell type assignment, stopping the algorithm.")
} else {
  print(paste0("Fraction of cells missing cell type assignment: ",round(mean(is.na(sample_metadata$celltype)),2)))
}


#######################
## Fetch ATAC matrix ##
#######################

if (file.exists(args$atac_peak_matrix)) {
  atac.se <- readRDS(args$atac_peak_matrix)[,sample_metadata$cell]
} else {

  # Load ArchR project
  source(here::here("atac/archR/load_archR_project.R"))
  stopifnot(sample_metadata$cell%in%rownames(ArchRProject))
  ArchRProject.filt <- ArchRProject[sample_metadata$cell]

  # Fetch ATAC peak matrix
  atac.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix", binarize = FALSE)
  dim(atac.se)

  # Define feature names
  row.ranges.dt <- rowRanges(atac.se) %>% as.data.table %>% 
    setnames("seqnames","chr") %>%
    .[,c("chr","start","end")] %>%
    .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
  rownames(atac.se) <- row.ranges.dt$idx
}

#######################
## Feature selection ##
#######################

# Load feature stats
atac_featureStats.dt <- fread(args$atac_feature_stats)
stopifnot(sort(atac_featureStats.dt$feature)==sort(rownames(atac.se)))

# Define highly variable features
atac_features <- atac_featureStats.dt %>% 
  setorder(-var_pseudobulk) %>% 
  head(n=opts$nfeatures) %>% .$feature

##############################
## ATAC TFIDF normalisation ##
##############################

stopifnot(atac_features%in%rownames(atac.se))

atac_tfidf.mtx <- tfidf(assay(atac.se[atac_features,]), method=1, scale.factor=1e4)

###########################
## Latent Semantic Index ##
###########################

svd <- irlba(atac_tfidf.mtx, opts$ndims, opts$ndims)
svdDiag <- matrix(0, nrow=opts$ndims, ncol=opts$ndims)
diag(svdDiag) <- svd$d
lsi.mtx <- t(svdDiag %*% t(svd$v))
rownames(lsi.mtx) <- colnames(atac.se)
colnames(lsi.mtx) <- paste0("LSI",seq_len(ncol(lsi.mtx)))

rm(svd,svdDiag); gc(reset=T)

###############
## Parse LSI ##
###############

# (Optional) Remove dimensions that correlate with nFrags
if (opts$remove_dim_cor_seq_depth) {
  lsi.mtx <- lsi.mtx[,which(cor(lsi.mtx,sample_metadata$nFrags_atac)<=0.75)]
}

# (Optional) scale dimensions, Z scores
# lsi.mtx <- sweep(lsi.mtx - rowMeans(lsi.mtx), 1, matrixStats::rowSds(lsi.mtx),`/`)
if (opts$scale_dims) {
  lsi.mtx <- sweep(lsi.mtx - colMeans(lsi.mtx), 2, matrixStats::colSds(lsi.mtx),`/`)
  lsi.mtx[lsi.mtx>2] <- 2
  lsi.mtx[lsi.mtx<(-2)] <- (-2)
}

##########
## UMAP ##
##########

umap.mtx <- uwot::umap(lsi.mtx, n_neighbors = 40, min_dist = 0.4, metric = "cosine", fast_sgd = TRUE) 

###########################################################################
## Create neighbourhood graph for cells that are missing cell type label ##
###########################################################################

cat("Creating neighbourhood graph...\n")

df.observed <- sample_metadata[!is.na(celltype),c("cell","celltype")]
cells.missing.label <- sample_metadata[is.na(celltype),cell]

X.observed <- lsi.mtx[df.observed$cell,]
X.missing <- lsi.mtx[cells.missing.label,]

# returns a list of 2 (N,K) matrices.:
# - nn.idx: 1-indexed indices
# - nn.dists: distances
knnObj <- nabor::knn(data = X.observed, query = X.missing, k = args$k)

####################
## kNN prediction ##
####################

cat("Doing kNN prediction...\n")

neighbour.celltypes <- apply(knnObj$nn.idx, 1, function(x) df.observed$celltype[x]) %>% t
predicted.celltype <- apply(neighbour.celltypes, 1, function(x) getmode(x, 1:length(x)))

############################
## Update sample metadata ##
############################

cat("Updating sample metadata...\n")

celltype_predictions.dt <- data.table(
  cell = cells.missing.label,
  celltype.predicted = predicted.celltype
)

tmp <- fread(args$metadata)
if (args$output_celltype_column%in%colnames(tmp)) tmp[[args$output_celltype_column]] <- NULL
sample_metadata.updated <- tmp %>% 
  merge(celltype_predictions.dt, by="cell", all.x=TRUE) %>%
  .[is.na(celltype.predicted),celltype.predicted:=eval(as.name(args$input_celltype_column))]

# parse metadata
sample_metadata.updated[pass_rnaQC==TRUE & pass_rnaQC==FALSE & !is.na(celltype.predicted),doublet_call:=FALSE]

# Save metadata
fwrite(sample_metadata.updated, file.path(args$outdir,"sample_metadata_after_celltype_assignment.txt.gz"), sep="\t", na="NA", quote=F)

#################
## Print stats ##
#################

print(sprintf("Number of cells that passed atac QC but did not have cell type assignment (before imputation): N=%s",tmp[pass_atacQC==TRUE & pass_rnaQC==FALSE,.N]))
print(sprintf("Number of cells that pass atac QC and now have cell type assignment (after imputation): N=%s",sample_metadata.updated[pass_atacQC==TRUE & pass_rnaQC==FALSE,.N]))

##########
## Plot ##
##########

cat("Plotting...\n")
    
to.plot <- umap.mtx %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("cell","umap1","umap2")) %>%
  merge(sample_metadata.updated, by="cell") %>%
  setnames(args$input_celltype_column,"celltype")

p1 <- ggplot(to.plot[!is.na(celltype)], aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype), size=1.25, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  labs(title=sprintf("Cells that have cell type assignment (%s, N=%s)",args$input_celltype_column,to.plot[!is.na(celltype),.N])) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=rel(0.9)))

p2 <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=is.na(celltype), size=is.na(celltype), alpha=is.na(celltype)), shape=21, color="black", stroke=0.05) +
  scale_size_manual(values=c("TRUE"=0.75, "FALSE"=0.4)) +
  scale_fill_manual(values=c("TRUE"="red", "FALSE"="gray60")) +
  scale_alpha_manual(values=c("TRUE"=0.75, "FALSE"=0.25)) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  labs(title=sprintf("Highlighting cells that do not have cell type assignment (in red, N=%s)",to.plot[is.na(celltype),.N])) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=rel(0.9)))
  
p3 <- ggplot(to.plot[is.na(celltype)], aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype.predicted), size=1.25, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  labs(title=sprintf("Subset of cells after celltype prediction (%s, N=%s)",args$output_celltype_column,to.plot[is.na(celltype),.N])) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=rel(0.9)))

p <- cowplot::plot_grid(plotlist=list(p1,p2,p3), nrow = 1)

pdf(file.path(args$outdir,"celltype_assignment.pdf"), width=16, height=6)
print(p)
dev.off()


