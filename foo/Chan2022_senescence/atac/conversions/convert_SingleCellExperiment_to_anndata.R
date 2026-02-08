suppressPackageStartupMessages(library(reticulate))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
io$sce <- file.path(io$basedir,"processed/atac/SingleCellExperiment.rds")
io$outfile <- paste0(io$basedir,"/processed/atac/anndata.h5ad")

#####################################
## Reticulate connection to scanpy ##
#####################################

use_condaenv("scanpy", required=TRUE)
sc <- import("scanpy")

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata)

###############
## Load data ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]

# Add sample metadata as colData
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

# Rename genes
# rownames(sce) <- rowData(sce)$gene_name

# Remove duplicated gene names
sce <- sce[!duplicated(rownames(sce))]

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

##########
## Save ##
##########

adata_sce$write_h5ad(io$outfile, compression="gzip")
