source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(UCell))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(BiocParallel))

set.seed(42)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outdir <- paste0(io$basedir,"/results/gene_signatures/cells"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$cores <- 8

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata)

table(cell_metadata.dt$assay)

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

##########################
## Load gene signatures ##
##########################

source(here::here("gene_signatures/load_gene_signatures.R"))

# tmp
tmp <- fread("/Users/rargelaguet/data/adult_fibroblast_reprogramming_10x_multiome_reik/results/rna/crazy_subpopulations/metacells/marker_genes_high_gapdh.txt.gz") %>% .[logFC<(-1.5)]
gene.sets <- list()
gene.sets[["gapdh_high"]] <- tmp$feature

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = T)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

###############
## Run UCell ##
###############

if (opts$cores>1) {
  BPPARAM <- BiocParallel::MulticoreParam(workers=opts$cores)
  sce <- ScoreSignatures_UCell(sce, features=gene.sets, BPPARAM=BPPARAM, maxRank=1000)
} else {
  sce <- ScoreSignatures_UCell(sce, features=gene.sets)  
}

# Save
saveRDS(altExp(sce,'UCell'), file.path(io$outdir,"UCell_SummarizedExperiment.rds"))
