suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))

#########
## I/O ##
#########

io <- list()

# Personal laptop
if (grepl("^(BAI-|SDI-|CAM-)",Sys.info()[['nodename']])) {
  if (Sys.info()[['user']]=="rargelaguet") {
    io$basedir <- "/Users/rargelaguet/data/Peyser2019_fibrosis"
    io$gene_metadata <- io$gene_metadata <- "/Users/rargelaguet/data/ensembl/mouse/v102/BioMart/all_genes/Mmusculus_genes_BioMart.102.txt"
    io$TFs <- "/Users/rargelaguet/data/mm10_regulation/TFs/TFs.txt"
    io$TEs <- "/Users/rargelaguet/data/hg38_regulation/transposable_elements/scTE/scTE_TEs.txt"
    io$cell_surface_markers <- "/Users/rargelaguet/data/hg38_regulation/surface_markers/cell_surface_markers.csv"
    io$epigenetic_factors <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/epigenetic_factors.txt"
    io$mito_genes <- "/Users/rargelaguet/data/mm10_regulation/gene_lists/mitochondrial_genes.txt"
    io$ribo_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/ribosome_genes.txt"
    io$dna_repair_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/dna_repair.txt"
    io$cell_cycle_g2m_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/cell_cycle_g2m.txt"
    io$histones_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/histones_genes.txt"
    io$ecm_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/ECM_genes.txt"
  } else {
    stop("User not recognised")
  }
  
# Dockyard
} else if (grepl("^(dockyard-)",Sys.info()[['nodename']])) { 
  io$basedir <- "/share/rargelaguet/data/Peyser2019_fibrosis"
  io$gene_metadata <- io$gene_metadata <- "/share/rargelaguet/data/ensembl/mouse/v102/BioMart/all_genes/Mmusculus_genes_BioMart.102.txt"
  io$TFs <- "/share/rargelaguet/data/mm10_regulation/TFs/TFs.txt"
  io$TEs <- "/share/rargelaguet/data/hg38_regulation/transposable_elements/scTE/scTE_TEs.txt"
  io$cell_surface_markers <- "/share/rargelaguet/data/hg38_regulation/surface_markers/cell_surface_markers.csv"
  io$epigenetic_factors <- "/share/rargelaguet/data/hg38_regulation/gene_lists/epigenetic_factors.txt"
  io$mito_genes <- "/share/rargelaguet/data/mm10_regulation/gene_lists/mitochondrial_genes.txt"
  io$ribo_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/ribosome_genes.txt"
  io$dna_repair_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/dna_repair.txt"
  io$cell_cycle_g2m_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/cell_cycle_g2m.txt"
  io$histones_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/histones_genes.txt"
  io$ecm_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/ECM_genes.txt"
} else {
  stop("Workstation not recognised")
}

io$metadata <- io$cell.metadata  <- io$cell_metadata <- file.path(io$basedir,"processed/cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$anndata <-file.path(io$basedir,"processed/anndata.h5ad")

#############
## Options ##
#############

opts <- list()

# opts$conditions <- c("Saline Control", "Bleomycin Fibrosis Model", "Bleomycin Fibrosis plus Nintedanib")
opts$conditions <- c("Saline_Control", "Bleomycin_day11", "Bleomycin_nintedanib_day11")

opts$celltypes <- c("dendritic_cells", "Type_II_alveolar", "Type_I_alveolar", "mesothelium", 
                    "fibroblasts", "endothelial_neutrophils", "macrophages", "endothelial", 
                    "proliferating_myeloid", "club_cells", "B_cells", "monocytes", 
                    "?", "proliferating_lymphocytes", "T_&_NK", "neutrophils")
