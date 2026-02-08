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
if (grepl("(BAI-.*|SDI-.*)",Sys.info()['nodename'])) {
  io$basedir <- "/Users/rargelaguet/data/Chan2022_senescence"
  io$gene_metadata <- "/Users/rargelaguet/data/ensembl/human/v110/Hsapiens_genes_BioMart.110.txt"
  io$TFs <- "/Users/rargelaguet/data/mm10_regulation/TFs/TFs.txt"
  io$cell_surface_markers <- "/Users/rargelaguet/data/hg38_regulation/surface_markers/cell_surface_markers.csv"
  io$epigenetic_factors <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/epigenetic_factors.txt"
  io$mito_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/mitochondrial_genes.txt"
  io$ribo_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/ribosome_genes.txt"
  io$dna_repair_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/dna_repair.txt"
  io$cell_cycle_g2m_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/cell_cycle_g2m.txt"
  io$histones_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/histones_genes.txt"
  io$ecm_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/ECM_genes.txt"
} else {
  io$basedir <- "/share/rargelaguet/data/Chan2022_senescence"
  io$gene_metadata <- "/home/rargelaguet/data/ensembl/human/v110/Hsapiens_genes_BioMart.110.txt"
  io$TFs <- "/share/rargelaguet/data/mm10_regulation/TFs/TFs.txt"
  io$cell_surface_markers <- "/share/rargelaguet/data/hg38_regulation/surface_markers/cell_surface_markers.csv"
  io$epigenetic_factors <- "/share/rargelaguet/data/hg38_regulation/gene_lists/epigenetic_factors.txt"
  io$mito_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/mitochondrial_genes.txt"
  io$ribo_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/ribosome_genes.txt"
  io$dna_repair_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/dna_repair.txt"
  io$cell_cycle_g2m_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/cell_cycle_g2m.txt"
  io$histones_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/histones_genes.txt"
  io$ecm_genes <- "/share/rargelaguet/data/hg38_regulation/gene_lists/ECM_genes.txt"
}

io$metadata <- io$sample_metadata <- file.path(io$basedir,"sample_metadata.txt")
io$sce <- io$se <- file.path(io$basedir,"rna/processed/SingleCellExperiment.rds")

#############
## Options ##
#############

opts <- list()

opts$samples <- c(
  "hTERT_TP1_rep1", "hTERT_TP1_rep2", "hTERT_TP1_rep3", 
  "hTERT_TP2_rep1", "hTERT_TP2_rep2", "hTERT_TP2_rep3", 
  "hTERT_TP4_rep1", "hTERT_TP4_rep2", "hTERT_TP4_rep3", 
  "hTERT_TP5_rep1", "hTERT_TP5_rep2", "hTERT_TP5_rep3", 
  "hTERT_TP6_rep1", "hTERT_TP6_rep2", "hTERT_TP6_rep3", 
  "hTERT_TP7_rep1", "hTERT_TP7_rep2", "hTERT_TP7_rep3", 
  "RS_PDL20_TP1_rep1", "RS_PDL20_TP1_rep2", "RS_PDL20_TP1_rep3", 
  "RS_PDL25_TP2_rep1", "RS_PDL25_TP2_rep2", "RS_PDL25_TP2_rep3", 
  "RS_PDL28_TP3_rep1", "RS_PDL28_TP3_rep2", "RS_PDL28_TP3_rep3", 
  "RS_PDL33_TP4_rep1", "RS_PDL33_TP4_rep2", "RS_PDL33_TP4_rep3", 
  "RS_PDL37_TP5_rep1", "RS_PDL37_TP5_rep2", "RS_PDL37_TP5_rep3", 
  "RS_PDL46_TP6_rep1", "RS_PDL46_TP6_rep2", "RS_PDL46_TP6_rep3", 
  "RS_PDL50_TP8_rep1", "RS_PDL50_TP8_rep2", "RS_PDL50_TP8_rep3", 
  # "RS_PDL45_TP7_rep1", "RS_PDL45_TP7_rep2", "RS_PDL45_TP7_rep3", 
  "RS_PDL52_TP9_rep1", "RS_PDL52_TP9_rep2", "RS_PDL52_TP9_rep3", 
  "RS_PDL53_TP10_rep1", "RS_PDL53_TP10_rep2", "RS_PDL53_TP10_rep3", 
  "CD_D10_rep1", "CD_D10_rep2", "CD_D10_rep3", 
  "CD_D1_rep1", "CD_D1_rep2", "CD_D1_rep3",
  "CD_D1.5_rep1", "CD_D1.5_rep2", "CD_D1.5_rep3",
  "CD_D2_rep1", "CD_D2_rep2", "CD_D2_rep3",
  "CD_D2.5_rep1", "CD_D2.5_rep2", "CD_D2.5_rep3",
  "CD_D2.75_rep1", "CD_D2.75_rep2", "CD_D2.75_rep3", 
  "CD_D3_rep1", "CD_D3_rep2", "CD_D3_rep3", 
  "CD_D3.5_rep1", "CD_D3.5_rep2", "CD_D3.5_rep3",
  "CD_D4_rep1", "CD_D4_rep2", "CD_D4_rep3",
  "CD_D7_rep1", "CD_D7_rep2", "CD_D7_rep3",
  "RIS_d2_no_xray_rep1", "RIS_d2_no_xray_rep2", "RIS_d2_no_xray_rep3", 
  "RIS_d3_xray_rep1", "RIS_d3_xray_rep2", "RIS_d3_xray_rep3", 
  "RIS_d4_xray_rep1", "RIS_d4_xray_rep2", "RIS_d4_xray_rep3", 
  "RIS_d5_xray_rep1", "RIS_d5_xray_rep2", "RIS_d5_xray_rep3", 
  "RIS_d6_xray_rep1", "RIS_d6_xray_rep2", "RIS_d6_xray_rep3", 
  "RIS_d9_xray_rep1", "RIS_d9_xray_rep2", "RIS_d9_xray_rep3"
)


opts$conditions <- c(
  "hTERT", 
  "Replicative_Senescence", 
  "Cell_Density", 
  "Radiation_Senescence"
)

opts$timepoints <- c(
  # "TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP7", "TP8", "TP9", "TP10", 
  "TP1", "TP2", "TP3", "TP4", "TP5", "TP6", "TP8", "TP9", "TP10", 
  "D1", "D1.5", "D2", "D2.5", "D2.75", "D3", "D3.5", "D4", "D7", "D5", "D6", "D9", "D10"
)
