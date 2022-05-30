here::i_am("rna/processing/extract_TFs_from_SingleCellExperiment.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',    type="character",    help='')
p$add_argument('--motif_annotation',    type="character",    help='')
# p$add_argument('--motif2gene',    type="character",    help='')
p$add_argument('--TF_file',    type="character",    help='')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# io$basedir <- file.path(io$basedir,"test")
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
# args$motif_annotation <- "JASPAR"
# # args$motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s_motif2gene.txt.gz",args$motif_annotation))
# args$TF_file <- "/Users/argelagr/data/mm10_regulation/TFs/TFs.txt"
# args$outfile <- file.path(io$basedir,sprintf("processed/rna/SingleCellExperiment_TFs_%s.rds",args$motif_annotation))
## END TEST

#########################
## Load RNA expression ##
#########################

rna.sce <- readRDS(args$sce)

################
## Subset TFs ##
################

# Load TF annotation
# motif2gene.dt <- fread(args$motif2gene) %>%
#   .[gene%in%toupper(rownames(rna.sce))]
# rna_tf.sce <- rna.sce[str_to_title(motif2gene.dt$gene),]

TFs <- fread(args$TF_file)[[1]]
TFs <- TFs[TFs%in%toupper(rownames(rna.sce))]

# Subset TFs
rna_tf.sce <- rna.sce[str_to_title(TFs),]
rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))

##########
## Save ##
##########

saveRDS(rna_tf.sce, args$outfile)
