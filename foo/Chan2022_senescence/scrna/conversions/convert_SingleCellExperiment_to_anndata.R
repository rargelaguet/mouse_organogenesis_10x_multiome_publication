suppressPackageStartupMessages({
  library("reticulate")
  library("SingleCellExperiment")
})

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$metadata <- file.path(io$basedir,"cell_metadata_standarised.txt.gz")
io$sce <- file.path(io$basedir,"processed/cellranger/SingleCellExperiment.rds")
io$outfile <- file.path(io$basedir,"processed/cellranger/anndata_G1.h5ad")

#####################################
## Reticulate connection to scanpy ##
#####################################

use_condaenv("scanpy", required=TRUE)
sc <- import("scanpy")

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[cell_cycle_phase%in%c("G1")]

table(cell_metadata.dt$cell_cycle_phase)
table(cell_metadata.dt$genetic_background)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- readRDS(io$sce)[,cell_metadata.dt$cell]

# update colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

# update rowData
rownames(sce) <- rowData(sce)$ens_id

#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata_sce <- sc$AnnData(
  X   = t(counts(sce)),
  obs = as.data.frame(colData(sce)),
  var = as.data.frame(rowData(sce)),
  # var = data.frame(gene=rownames(sce), transcript_type=rowData(sce)$transcript_type, row.names=rownames(sce)),
  dtype = 'float32'
)

adata_sce

##########
## Save ##
##########

adata_sce$write_h5ad(io$outfile, compression="gzip")
