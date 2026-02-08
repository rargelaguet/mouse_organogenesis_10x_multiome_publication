suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))

#########
## I/O ##
#########

io <- list()
if (grepl("(BAI-.*|SDI-.*)",Sys.info()['nodename'])) {
  io$basedir <- "/Users/rargelaguet/data/Chan2022_senescence"
  io$gene_metadata <- "/Users/rargelaguet/data/ensembl/human/v99/bioMart/all_genes/Hsapiens_genes_BioMart.99.txt.gz"
} else {
  io$basedir <- "/share/rargelaguet/data/Chan2022_senescence"
  io$gene_metadata <- "/share/rargelaguet/data/ensembl/human/v99/bioMart/all_genes/Hsapiens_genes_BioMart.99.txt.gz"
}

io$metadata <- io$sample_metadata <- paste0(io$basedir,"/processed/atac/sample_metadata.txt")
io$atac.se <- file.path(io$basedir,"processed/gene_SummarizedExperiment.rds")

#############
## Options ##
#############

opts <- list()

opts$samples <- c(
  "pA_Control_REP1",
  "pA_Control_REP2",
  "pA_Control_REP3",
  "pA_TG_REP1",
  "pA_TG_REP2",
  "pA_TG_REP3",
  "pB_Control_REP1",
  "pB_Control_REP2",
  "pB_Control_REP3",
  "pB_TG_REP1",
  "pB_TG_REP2",
  "pB_TG_REP3",
  "pC_Control_REP1",
  "pC_Control_REP2",
  "pC_Control_REP3",
  "pC_TG_REP1",
  "pC_TG_REP2",
  "pC_TG_REP3",
  "pD_Control_REP1",
  "pD_Control_REP2",
  "pD_Control_REP3",
  "pD_TG_REP1",
  "pD_TG_REP2",
  "pD_TG_REP3",
  "pE_Control_REP1",
  "pE_Control_REP2",
  "pE_TG_REP1",
  "pE_TG_REP2",
  "pF_Control_REP1",
  "pF_Control_REP2",
  "pF_Control_REP3",
  "pF_TG_REP1",
  "pF_TG_REP2",
  "pF_TG_REP3",
  "pG_Control_REP1",
  "pG_Control_REP2",
  "pG_Control_REP3",
  "pG_TG_REP1",
  "pG_TG_REP2"
)

opts$sample2alias <- c(
  "pA_Control_REP1" = "Control_p8_rep1",
  "pA_Control_REP2" = "Control_p8_rep2",
  "pA_Control_REP3" = "Control_p8_rep3",
  "pA_TG_REP1" = "TG_p8_rep1",
  "pA_TG_REP2" = "TG_p8_rep2",
  "pA_TG_REP3" = "TG_p8_rep3",
  "pB_Control_REP1" = "Control_p11_rep1",
  "pB_Control_REP2" = "Control_p11_rep2",
  "pB_Control_REP3" = "Control_p11_rep3",
  "pB_TG_REP1" = "TG_p11_rep1",
  "pB_TG_REP2" = "TG_p11_rep2",
  "pB_TG_REP3" = "TG_p11_rep3",
  "pC_Control_REP1" = "p14_Control_rep1",
  "pC_Control_REP2" = "p14_Control_rep2",
  "pC_Control_REP3" = "p14_Control_rep3",
  "pC_TG_REP1" = "p14_TG_rep1",
  "pC_TG_REP2" = "p14_TG_rep2",
  "pC_TG_REP3" = "p14_TG_rep3",
  "pD_Control_REP1" = "Control_p19_rep1",
  "pD_Control_REP2" = "Control_p19_rep2",
  "pD_Control_REP3" = "Control_p19_rep3",
  "pD_TG_REP1" = "TG_p19_rep1",
  "pD_TG_REP2" = "TG_p19_rep2",
  "pD_TG_REP3" = "TG_p19_rep3",
  "pE_Control_REP1" = "Control_p22_rep1",
  "pE_Control_REP2" = "Control_p22_rep2",
  "pE_TG_REP1" = "TG_p22_rep1",
  "pE_TG_REP2" = "TG_p22_rep2",
  "pF_Control_REP1" = "Control_p24_rep1",
  "pF_Control_REP2" = "Control_p24_rep2",
  "pF_Control_REP3" = "Control_p24_rep3",
  "pF_TG_REP1" = "TG_p24_rep1",
  "pF_TG_REP2" = "TG_p24_rep2",
  "pF_TG_REP3" = "TG_p24_rep3",
  "pG_Control_REP1" = "Control_p29_rep1",
  "pG_Control_REP2" = "Control_p29_rep2",
  "pG_Control_REP3" = "Control_p29_rep3",
  "pG_TG_REP1" = "TG_p29_rep1",
  "pG_TG_REP2" = "TG_p29_rep2",
)

# to filter
opts$aliases <- c(
  "Control_p8_rep1",
  "Control_p8_rep2",
  "Control_p8_rep3",
  "Control_p11_rep1",
  "Control_p11_rep2",
  "Control_p11_rep3",
  "Control_p14_rep1",
  "Control_p14_rep2",
  "Control_p14_rep3",
  "Control_p19_rep1",
  "Control_p19_rep2",
  "Control_p19_rep3",
  "Control_p22_rep1",
  "Control_p22_rep2",
  "Control_p22_rep3",
  "Control_p24_rep1",
  "Control_p24_rep2",
  "Control_p24_rep3",
  "Control_p29_rep1",
  "Control_p29_rep2",
  "Control_p29_rep3",
  "TG_p8_rep1",
  "TG_p8_rep2",
  "TG_p8_rep3",
  "TG_p11_rep1",
  "TG_p11_rep2",
  "TG_p11_rep3",
  "TG_p14_rep1",
  "TG_p14_rep2",
  "TG_p14_rep3",
  "TG_p19_rep1",
  "TG_p19_rep2",
  "TG_p19_rep3",
  "TG_p22_rep1",
  "TG_p22_rep2",
  "TG_p22_rep3",
  "TG_p24_rep1",
  "TG_p24_rep2",
  "TG_p24_rep3",
  "TG_p29_rep1",
  "TG_p29_rep2",
  "TG_p29_rep3"
)

opts$passages <- c(8L, 11L, 14L, 19L, 22L, 24L, 29L)

opts$treatments <- c("Control", "TG")

opts$groups <- c(
  "p8_Control", 
  "p11_Control", 
  "p14_Control", 
  "p19_Control", 
  "p22_Control", 
  "p24_Control", 
  "p29_Control", 
  "p8_TG", 
  "p11_TG", 
  "p14_TG", 
  "p19_TG", 
  "p22_TG", 
  "p24_TG", 
  "p29_TG"
)
