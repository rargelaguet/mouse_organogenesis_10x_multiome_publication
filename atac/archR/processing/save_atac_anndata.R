here::i_am("atac/archR/processing/save_atac_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(reticulate))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--python',    type="character",    help='')
p$add_argument('--atac_matrix',    type="character",    help='')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--outfile',     type="character",  help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$python = "/Users/argelagr/opt/anaconda3/envs/main/bin/python" # "/bi/group/reik/ricard/software/miniconda3/envs/main/bin/python"
# args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$atac_matrix <- file.path(io$basedir,"processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
# args$outfile <- file.path(io$basedir,"processed/atac/anndata/PeakMatrox_anndata.h5ad")
## END TEST ##

dir.create(dirname(args$outfile), showWarnings = F, recursive = T)

#####################################
## Reticulate connection to scanpy ##
#####################################

use_python(args$python, required=TRUE)
sc <- import("scanpy")

#####################
## Define settings ##
#####################

##########################
## Load sample metadata ##
##########################

metadata.dt <- fread(args$metadata) %>% 
  .[pass_rnaQC==TRUE & pass_atacQC==TRUE & doublet_call==FALSE & !is.na(celltype)] %>%
  .[,c("cell","sample","stage","genotype","celltype","nFrags_atac","nFeature_RNA")]# %>%
  # setnames("celltype.predicted","celltype")

fwrite(metadata.dt, file.path(dirname(args$outfile),"cell_metadata.txt.gz"), sep="\t", quote=F, na="NA")

######################
## Load atac matrix ##
######################

atac.se <- readRDS(args$atac_matrix)[,metadata.dt$cell]

#############################################
## Convert SingleCellExperiment to AnnData ##
#############################################

adata <- sc$AnnData(
    X   = t(assay(atac.se)),
    obs = as.data.frame(colData(atac.se)),
    var = as.data.frame(rowData(atac.se))
)
print(adata)
print(head(adata$obs))
print(head(adata$var))

##########################
## Parse anndata object ##
##########################

adata$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata$obs$celltype)))])
adata$uns$update(stage_colors = opts$stage.colors[sort(unique(as.character(adata$obs$stage)))])

##########
## Save ##
##########

adata$write_h5ad(args$outfile)
