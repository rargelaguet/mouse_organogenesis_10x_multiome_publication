source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(UCell))
suppressPackageStartupMessages(library(msigdbr))

set.seed(42)

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- file.path(io$basedir,"rna/processed/sample_metadata.txt")
# io$sce <- file.path(io$basedir,"rna/processed/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"rna/results/gene_signatures/UCell"); dir.create(io$outdir, showWarnings = F, recursive = T)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] %>% .[,condition:=factor(condition,levels=opts$conditions)] %>%
  .[time%in%opts$timepoints] %>% .[,time:=factor(time,levels=opts$timepoints)] %>% droplevels

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

###############
## Load data ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]
rownames(sce) <- rowData(sce)$gene_name
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

##########################
## Load gene signatures ##
##########################

source(here::here("rna/gene_signatures/load_gene_signatures.R"))

# myofibroblast_signature_champ.dt <- fread("/Users/rargelaguet/data/hg38_regulation/CHAMP_signatures/champMyoSig.csv", header = T, select=c(2,3)) %>%
#   setnames(c("gene","ens_id"))
# 
# gapdh_signature.dt <- fread("/Users/rargelaguet/data/fibros_reik_lab/results/rna/crazy_subpopulations/metacells/marker_genes_high_gapdh.txt.gz", select=c(1,2,6)) %>%
#   .[transcript_type=="gene"] %>% .[,transcript_type:=NULL]  %>% .[logFC<(-2)]
# 
# gene.sets[["myofibroblast_champ"]] <- myofibroblast_signature_champ.dt$gene
# gene.sets[["gapdh_signature"]] <- gapdh_signature.dt$feature

# yamanaka_regulons <- readRDS("/Users/rargelaguet/data/Rand2021_oskm_combinations_microarray/results/regulons/yamanaka_regulons_100_Rand2017.rds")
# gene.sets <- yamanaka_regulons
# names(gene.sets) <- paste0(names(gene.sets),"_regulon")

# grep("senescence",names(gene.sets), ignore.case = TRUE, value=T)
gene.sets <- gene.sets[grep("senescence",names(gene.sets), ignore.case = TRUE)]

names(gene.sets)

###############
## Run UCell ##
###############

sce <- ScoreSignatures_UCell(sce, features=gene.sets, maxRank=1e4)

# Save
# saveRDS(altExp(sce,'UCell'), file.path(io$outdir,"UCell_SingleCellExperiment.rds"))
saveRDS(altExp(sce,'UCell'), file.path(io$outdir,"UCell_SingleCellExperiment_senescence.rds"))
# saveRDS(altExp(sce,'UCell'), file.path(io$outdir,"UCell_SingleCellExperiment_Yamanaka_regulons.rds"))
