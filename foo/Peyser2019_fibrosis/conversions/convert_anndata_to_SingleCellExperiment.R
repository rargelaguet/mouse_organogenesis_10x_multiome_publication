source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(scuttle))
suppressPackageStartupMessages(library(zellkonverter))
suppressPackageStartupMessages(library(Matrix))

#####################
## Define settings ##
#####################

io$python_path <- "/Users/rargelaguet/miniforge3/envs/scanpy/bin/python"
# io$python_path <- "/share/rargelaguet/mambaforge/envs/scanpy/bin/python"
io$anndata <- file.path(io$basedir,"processed/anndata.h5ad")
io$outdir <- file.path(io$basedir,"processed")

################
## Reticulate ##
################

reticulate::use_python(io$python_path, required = TRUE)

sc <- import("scanpy")

############################################
## Load anndata into SingleCellExperiment ##
############################################

# sce <- readH5AD(io$anndata, X_name="raw", use_hdf5 = FALSE, reader = "python")
sce <- readH5AD(io$anndata, X_name="raw", use_hdf5 = FALSE, reader = "R")

print("Overview of colData")
head(colData(sce))

print("Overview of rowData")
head(rowData(sce))

# Set gene names
if (is.null(rownames(sce))) {
  if ("gene"%in%colnames(rowData(sce))) {
    rownames(sce) <- rowData(sce)$gene
  }
}
stopifnot(!is.null(rownames(sce)))
print("Overview of gene names")
head(rownames(sce))

# Set cell names
if (is.null(colnames(sce))) {
  if ("cell"%in%colnames(colData(sce))) {
    colnames(sce) <- colData(sce)$cell
  }
}
stopifnot(!is.null(colnames(sce)))
print("Overview of cell names")
head(colnames(sce))

# Assays
print("Assays:")
print(assayNames(sce))

# Subset assays
assays(sce) <- assays(sce)["raw"]

# set assay names
assayNames(sce) <- "counts"
# assayNames(sce) <- c("counts","spliced","unspliced")

print("Overview of counts")
counts(sce)[1:10,1:10]

# Convert counts to sparse matrix if not
if (!is(counts(sce), 'sparseMatrix')) {
  counts(sce) <- Matrix(counts(sce), sparse=TRUE)
}

# reducedDims
print(reducedDims(sce))

# extract cell metadata
cell_metadata.dt <- colData(sce) %>% as.data.table(keep.rownames=T) %>%
  setnames("rn","cell")

###########
## Parse ##
###########

# fwrite(cell_metadata.dt, file.path(io$outdir,"cell_metadata.txt.gz"), sep="\t", na="NA", quote=F)
saveRDS(sce, file.path(io$outdir,"SingleCellExperiment.rds"))

