suppressPackageStartupMessages(library(reticulate))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$metadata <- file.path(io$basedir,"rna/sample_metadata_standarised.txt")
io$sce <- file.path(io$basedir,"rna/processed/SingleCellExperiment.rds")
io$outfile <- file.path(io$basedir,"rna/processed/anndata_standarised.h5ad")

#####################################
## Reticulate connection to scanpy ##
#####################################

use_condaenv("scanpy", required=TRUE)
sc <- import("scanpy")

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata)

colnames(sample_metadata.dt) <- colnames(sample_metadata.dt) %>% gsub(" ","_",.)

###############
## Load data ##
###############

# sce <- readRDS(io$sce)[,sample_metadata.dt$sample]
sce <- readRDS(io$sce)[,sample_metadata.dt$sample_name]

# Add sample metadata as colData
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample_name") %>% DataFrame

# Rename genes
rowData(sce)$ens_id <- rownames(sce)
# rownames(sce) <- rowData(sce)$gene_name

# Remove duplicated gene names
sce <- sce[!duplicated(rownames(sce))]

###############
## Subset data ##
###############

sce <- sce[,sce$other_perturbation%in%c("hTERT", "Replicative_Senescence", "Radiation_Senescence")]

#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata_sce <- sc$AnnData(
    X   = t(counts(sce)),
    obs = as.data.frame(colData(sce)),
    var = data.frame(gene=rownames(sce), ens_id=rowData(sce)$ens_id, row.names=rownames(sce)),
    dtype=  'float32'
)

adata_sce

##########
## Save ##
##########

adata_sce$write_h5ad(io$outfile, compression="gzip")
