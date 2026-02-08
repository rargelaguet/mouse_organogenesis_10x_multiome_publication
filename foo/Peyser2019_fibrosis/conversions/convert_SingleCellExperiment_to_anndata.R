suppressPackageStartupMessages({
  library("reticulate")
  library("SingleCellExperiment")
})

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# cellranger
io$metadata <- file.path(io$basedir,"processed/cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outfile <- file.path(io$basedir,"processed/anndata.h5ad")

#####################################
## Reticulate connection to scanpy ##
#####################################

use_condaenv("scanpy", required=TRUE)
sc <- import("scanpy")

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

# table(sce$sample)
  
#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata_sce <- sc$AnnData(
    X   = t(counts(sce)),
    obs = as.data.frame(colData(sce)),
    var = data.frame(gene=rownames(sce), row.names=rownames(sce)),
    dtype=  'float32'
)

adata_sce

# Add precomputed dimensionality reduction
adata_sce$obsm$update(umap = reducedDim(sce, "UMAP"))
adata_sce$obsm$update(pca = reducedDim(sce, "PCA"))

##########
## Save ##
##########

adata_sce$write_h5ad(io$outfile, compression="gzip")
