here::i_am("rna/scanpy/create_anndata_from_SingleCellExperiment.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Load libraries
suppressPackageStartupMessages({
  library("reticulate")
  library("SingleCellExperiment")
})

#####################
## Define settings ##
#####################

io$outfile <- file.path(io$basedir,"processed/rna/anndata.h5ad")

#####################################
## Reticulate connection to scanpy ##
#####################################

sc <- import("scanpy")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype.mapped)]

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame


#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata_sce <- sc$AnnData(
    X   = t(counts(sce)),
    obs = as.data.frame(colData(sce)),
    var = data.frame(gene=rownames(sce), row.names=rownames(sce))
)
# adata_sce$obsm$update(umap = reducedDim(sce, "umap"))

adata_sce

# Add cell type colors
# colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['celltype']))]
# adata.uns['celltype'] = colPalette_celltypes
# colPalette_stages = [opts["stage_colors"][i.replace(" ","_").replace("/","_")] for i in sorted(np.unique(adata.obs['stage']))]
# adata.uns['stage_colors'] = colPalette_stages
adata_sce$uns$update(celltype.mapped_colors = opts$celltype.colors[sort(unique(as.character(adata_sce$obs$celltype.mapped)))])
adata_sce$uns$update(stage_colors = opts$stage.colors[sort(unique(as.character(adata_sce$obs$stage)))])
adata_sce$uns["celltype.mapped_colors"]
adata_sce$uns["stage_colors"]

##########
## Save ##
##########

adata_sce$write_h5ad(io$outfile)
