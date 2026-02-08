suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(Seurat))

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sample',             type="character")

args <- p$parse_args(commandArgs(TRUE))

# Test
# args$sample <- c("SRR11427626_GSM4433529_D0")

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outdir <- paste0(io$basedir,"/results/cell_cycle"); dir.create(io$outdir, showWarnings = F, recursive=TRUE)

# Sanity checks
# stopifnot(args$sample%in%opts$samples)
# stopifnot(length(args$sample)==1)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata)# %>% .[sample%in%args$sample]

table(cell_metadata.dt$sample)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

############
## Seurat ##
############

seurat <- as.Seurat(sce)

# Normalisation
# seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

# Feature selection
# seurat <- FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = 2500)

# Scaling
seurat <- ScaleData(seurat)

###########################
## Cell cycle assignment ##
###########################

s.genes = gprofiler2::gorth(cc.genes$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
g2m.genes = gprofiler2::gorth(cc.genes$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes

seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

cell_cycle.dt <- seurat@meta.data[,c("S.Score","G2M.Score","Phase")] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("cell","cell_cycle_s_score","cell_cycle_g2m_score","cell_cycle_phase"))

cell_cycle.dt[,cell_cycle_s_score:=round(cell_cycle_s_score,2)]
cell_cycle.dt[,cell_cycle_g2m_score:=round(cell_cycle_g2m_score,2)]

##########
## Plot ##
##########

# RidgePlot(seurat, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# 
# seurat <- RunPCA(seurat, features = c(s.genes, g2m.genes))
# DimPlot(seurat)

##########################
## Update cell metadata ##
##########################

cell_metadata_cell_cycle.dt <- cell_metadata.dt %>% 
  merge(cell_cycle.dt, by=c("cell"), all.x=TRUE)

##########
## Save ##
##########

# fwrite(cell_cycle.dt, file.path(io$outdir,sprintf("%s_cell_cycle.txt.gz",args$sample)), sep="\t", quote=F)
# fwrite(cell_cycle.dt, file.path(io$outdir,"cell_cycle_assignments.txt.gz"), sep="\t", quote=F)
fwrite(cell_metadata_cell_cycle.dt, file.path(io$outdir,"cell_metadata_after_cell_cycle.txt.gz"), sep="\t", quote=F, na="NA")



