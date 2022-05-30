here::i_am("rna/conversions/convert_SingleCellExperiment_to_anndata.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(scuttle))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--python_path',   type="character",    help='Python path for reticulate')
p$add_argument('--metadata',   type="character",    help='Cell metadata')
p$add_argument('--sce',  type="character",              help='SingleCellExperiment input file') 
p$add_argument('--outfile',          type="character",                help='Anndata output file')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset data')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$python_path <- "/bi/group/reik/ricard/software/miniconda3/envs/main/bin/python" # "/Users/argelagr/opt/anaconda3/envs/main/bin/python"
# # args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
# args$outfile <- file.path(io$basedir,"processed/rna/anndata.h5ad")
# args$test_mode <- FALSE
## END TEST ##

################
## Reticulate ##
################

reticulate::use_python(args$python_path, required = TRUE)

sc <- import("scanpy")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE]# %>%
  # .[,c("cell", "sample", "stage", "nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA", "celltype.mapped")] %>%
  # setnames("celltype.mapped","celltype")

if (args$test_mode) {
	print("Test mode activated, subsetting number of cells...")
	sample_metadata <- sample_metadata %>% head(n=100)
}

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

###########
## Parse ##
###########

# reducedDims
# reducedDimNames(sce)[1] <- "X_pca_precomputed"
# forceatlas.mtx <- colData(sce)[,c("cell","Dim1","Dim2")] %>% matrix.please
# reducedDims(sce)[["X_draw_graph_fa_precomputed"]] <- forceatlas.mtx

# colData
# colData(sce) <- colData(sce)[,c("cell","Dim1","Dim2","sample", "stage.mapped", "somite_count", "tube_name", "tube_name_corrected", "celltype.clustering")]
# colnames(colData(sce)) <- c("cell","x","y","sample", "stage", "somite_count", "tube_name", "tube_name_corrected", "celltype")


#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

# adata_sce <- sc$AnnData(
#     X   = Matrix::t(counts(sce)),
#     obs = as.data.frame(colData(sce)),
#     obsm = list(
#       "X_pca_precomputed" = as.matrix(reducedDim(sce, "X_pca_precomputed")),
#       "X_draw_graph_fa_precomputed" = as.matrix(reducedDim(sce, "X_draw_graph_fa_precomputed"))
#     ),
#     var = data.frame(gene=rownames(sce), row.names=rownames(sce))
# )

adata_sce <- sc$AnnData(
    X   = t(counts(sce)),
    obs = as.data.frame(colData(sce)),
    var = data.frame(gene=rownames(sce), row.names=rownames(sce))
)

adata_sce

# TO-DO: CAN WE MAKE THE MATRICES SPARSE??

##########################
## Parse anndata object ##
##########################

# Add stage colors
adata_sce$uns$update(stage_colors = opts$stage.colors[sort(unique(as.character(adata_sce$obs$stage)))])
adata_sce$uns["stage_colors"]

# Add celltype colors
adata_sce$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata_sce$obs$celltype)))])
adata_sce$uns["celltype_colors"]

##########
## Save ##
##########

adata_sce$write_h5ad(args$outfile, compression="gzip")
