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
if (grepl("^(BAI-|SDI-|CAM-)", Sys.info()['nodename'])) {
  io$basedir <- "/Users/rargelaguet/data/Chan2022_senescence/scrna"
  io$gene_metadata <- "/Users/rargelaguet/data/ensembl/human/v99/bioMart/all_genes/Hsapiens_genes_BioMart.99.txt.gz"
} else {
  io$basedir <- "/home/rargelaguet/data/Chan2022_senescence/scrna"
  io$gene_metadata <- "/home/rargelaguet/data/ensembl/human/v99/bioMart/all_genes/Hsapiens_genes_BioMart.99.txt.gz"
}
io$metadata <- io$cell_metadata <- file.path(io$basedir,"cell_metadata.txt.gz")
io$sample_metadata <- file.path(io$basedir,"sample_metadata.txt")
io$sce <- io$rna.sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$anndata <- file.path(io$basedir,"processed/anndata.h5ad")

#############
## Options ##
#############

opts <- list()

opts$samples <- c( 
  "sample-SRR14646407_GSM5341090_RS_sc_hTert_t2",
  "sample-SRR14646408_GSM5341091_RS_sc_hTert_t4",
  "sample-SRR14646409_GSM5341092_RS_sc_hTert_t5",
  "sample-SRR14646410_GSM5341093_RS_sc_hTert_t6",
  "sample-SRR14646411_GSM5341094_RS_sc_hTert_t7",
  "sample-SRR14646412_GSM5341095_RS_sc_PDL_25_1",
  "sample-SRR14646413_GSM5341096_RS_sc_PDL_29_1",
  "sample-SRR14646414_GSM5341097_RS_sc_PDL_29_2",
  "sample-SRR14646415_GSM5341098_RS_sc_PDL_33_1",
  "sample-SRR14646416_GSM5341099_RS_sc_PDL_33_2",
  "sample-SRR14646417_GSM5341100_RS_sc_PDL_37_1",
  "sample-SRR14646418_GSM5341101_RS_sc_PDL_37_2",
  "sample-SRR14646419_GSM5341102_RS_sc_PDL_46_1",
  "sample-SRR14646420_GSM5341103_RS_sc_PDL_46_2",
  "sample-SRR14646421_GSM5341104_RS_sc_PDL_50_1",
  "sample-SRR14646422_GSM5341105_RS_sc_PDL_50_2"
)

opts$aliases <- c( 
  "hTert_t2",
  "hTert_t4",
  "hTert_t5",
  "hTert_t6",
  "hTert_t7",
  "PDL_25_rep1",
  "PDL_29_rep1",
  "PDL_29_rep2",
  "PDL_33_rep1",
  "PDL_33_rep2",
  "PDL_37_rep1",
  "PDL_37_rep2",
  "PDL_46_rep1",
  "PDL_46_rep2",
  "PDL_50_rep1",
  "PDL_50_rep2"
)

opts$groups <- c(
  "hTert_t2", "hTert_t4", "hTert_t5", "hTert_t6", "hTert_t7", 
  "PDL_25", "PDL_29", "PDL_33", "PDL_37", "PDL_46", "PDL_50"
)

opts$sample2alias <- c( 
  "sample-SRR14646407_GSM5341090_RS_sc_hTert_t2" =  "hTert_t2",
  "sample-SRR14646408_GSM5341091_RS_sc_hTert_t4" =  "hTert_t4",
  "sample-SRR14646409_GSM5341092_RS_sc_hTert_t5" =  "hTert_t5",
  "sample-SRR14646410_GSM5341093_RS_sc_hTert_t6" =  "hTert_t6",
  "sample-SRR14646411_GSM5341094_RS_sc_hTert_t7" =  "hTert_t7",
  "sample-SRR14646412_GSM5341095_RS_sc_PDL_25_1" =  "PDL_25_rep1",
  "sample-SRR14646413_GSM5341096_RS_sc_PDL_29_1" =  "PDL_29_rep1",
  "sample-SRR14646414_GSM5341097_RS_sc_PDL_29_2" =  "PDL_29_rep2",
  "sample-SRR14646415_GSM5341098_RS_sc_PDL_33_1" =  "PDL_33_rep1",
  "sample-SRR14646416_GSM5341099_RS_sc_PDL_33_2" =  "PDL_33_rep2",
  "sample-SRR14646417_GSM5341100_RS_sc_PDL_37_1" =  "PDL_37_rep1",
  "sample-SRR14646418_GSM5341101_RS_sc_PDL_37_2" =  "PDL_37_rep2",
  "sample-SRR14646419_GSM5341102_RS_sc_PDL_46_1" =  "PDL_46_rep1",
  "sample-SRR14646420_GSM5341103_RS_sc_PDL_46_2" =  "PDL_46_rep2",
  "sample-SRR14646421_GSM5341104_RS_sc_PDL_50_1" =  "PDL_50_rep1",
  "sample-SRR14646422_GSM5341105_RS_sc_PDL_50_2" =  "PDL_50_rep2"
)


opts$sample2group <- c( 
  "sample-SRR14646407_GSM5341090_RS_sc_hTert_t2" =  "hTert_t2",
  "sample-SRR14646408_GSM5341091_RS_sc_hTert_t4" =  "hTert_t4",
  "sample-SRR14646409_GSM5341092_RS_sc_hTert_t5" =  "hTert_t5",
  "sample-SRR14646410_GSM5341093_RS_sc_hTert_t6" =  "hTert_t6",
  "sample-SRR14646411_GSM5341094_RS_sc_hTert_t7" =  "hTert_t7",
  "sample-SRR14646412_GSM5341095_RS_sc_PDL_25_1" =  "PDL_25",
  "sample-SRR14646413_GSM5341096_RS_sc_PDL_29_1" =  "PDL_29",
  "sample-SRR14646414_GSM5341097_RS_sc_PDL_29_2" =  "PDL_29",
  "sample-SRR14646415_GSM5341098_RS_sc_PDL_33_1" =  "PDL_33",
  "sample-SRR14646416_GSM5341099_RS_sc_PDL_33_2" =  "PDL_33",
  "sample-SRR14646417_GSM5341100_RS_sc_PDL_37_1" =  "PDL_37",
  "sample-SRR14646418_GSM5341101_RS_sc_PDL_37_2" =  "PDL_37",
  "sample-SRR14646419_GSM5341102_RS_sc_PDL_46_1" =  "PDL_46",
  "sample-SRR14646420_GSM5341103_RS_sc_PDL_46_2" =  "PDL_46",
  "sample-SRR14646421_GSM5341104_RS_sc_PDL_50_1" =  "PDL_50",
  "sample-SRR14646422_GSM5341105_RS_sc_PDL_50_2" =  "PDL_50"
)

opts$sample2day <- c( 
  "sample-SRR14646407_GSM5341090_RS_sc_hTert_t2" =  "25",
  "sample-SRR14646408_GSM5341091_RS_sc_hTert_t4" =  "29",
  "sample-SRR14646409_GSM5341092_RS_sc_hTert_t5" =  "33",
  "sample-SRR14646410_GSM5341093_RS_sc_hTert_t6" =  "37",
  "sample-SRR14646411_GSM5341094_RS_sc_hTert_t7" =  "46",
  "sample-SRR14646412_GSM5341095_RS_sc_PDL_25_1" =  "25",
  "sample-SRR14646413_GSM5341096_RS_sc_PDL_29_1" =  "29",
  "sample-SRR14646414_GSM5341097_RS_sc_PDL_29_2" =  "29",
  "sample-SRR14646415_GSM5341098_RS_sc_PDL_33_1" =  "33",
  "sample-SRR14646416_GSM5341099_RS_sc_PDL_33_2" =  "33",
  "sample-SRR14646417_GSM5341100_RS_sc_PDL_37_1" =  "37",
  "sample-SRR14646418_GSM5341101_RS_sc_PDL_37_2" =  "37",
  "sample-SRR14646419_GSM5341102_RS_sc_PDL_46_1" =  "46",
  "sample-SRR14646420_GSM5341103_RS_sc_PDL_46_2" =  "46",
  "sample-SRR14646421_GSM5341104_RS_sc_PDL_50_1" =  "50",
  "sample-SRR14646422_GSM5341105_RS_sc_PDL_50_2" =  "50"
)