here::i_am("rna/conversions/convert_SingleCellExperiment_to_anndata.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(scuttle))
library(zellkonverter)

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--python_path',   type="character",    help='Python path for reticulate')
p$add_argument('--anndata',  type="character",              help='Anndata input file') 
p$add_argument('--outfile',          type="character",                help='SingleCellExperiment output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args <- list()
args$python_path <- "/bi/group/reik/ricard/software/miniconda3/envs/main/bin/python" # "/Users/argelagr/opt/anaconda3/envs/main/bin/python"
args$anndata <- file.path(io$basedir,"processed/rna/velocyto/anndata_velocyto.h5ad")
args$outfile <- file.path(io$basedir,"processed/rna/velocyto/SingleCellExperiment_velocyto.rds")
## END TEST ##

################
## Reticulate ##
################

reticulate::use_python(args$python_path, required = TRUE)

sc <- import("scanpy")

############################################
## Load anndata into SingleCellExperiment ##
############################################

sce <- readH5AD(args$anndata, use_hdf5 = FALSE, reader = "python")

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

# set assay names
# assayNames(sce) <- "counts"
assayNames(sce) <- c("counts","spliced","unspliced")
print("Overview of counts")
counts(sce)[1:10,1:10]

# reducedDims

###########
## Parse ##
###########

saveRDS(sce, args$outfile)
