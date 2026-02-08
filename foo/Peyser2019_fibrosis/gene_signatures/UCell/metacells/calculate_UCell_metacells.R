source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(UCell))

set.seed(42)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"results/mouse/metacells/per_tissue/metacells_metadata.txt.gz")
io$sce <- file.path(io$basedir,"results/mouse/metacells/per_tissue/SingleCellExperiment_metacells.rds")
io$outdir <- paste0(io$basedir,"/results/mouse/gene_signatures/metacells"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$tissues <- c(
  "Bone", "Visceral Adipose", "Omentum", "Intestine", "Heart", 
  "Subcutaneous Adipose", "Liver", "Lymph Node", "Lung", "Mesentery", 
  "Muscle", "Pancreas", "Skin", "Spleen", "Tendon", "Artery"
)

opts$clusters <- c(
  "Ccl19", "Cxcl12", "Comp", "Col15a1", "Coch", "Pi16", "Fbln1", "Npnt", "Hhip", "Bmp4"
)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[Tissue%in%opts$tissues] %>%
  .[ClustName%in%opts$clusters] %>%
  setnames("metacell","cell")

table(cell_metadata.dt$ClustName)
table(cell_metadata.dt$Tissue)

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata_mouse)

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = F)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

##########################
## Load gene signatures ##
##########################

source(here::here("mouse/gene_signatures/load_gene_signatures.R"))

# gene.sets <- gene.sets[grep("grn",names(gene.sets))]
# gene.sets <- sapply(gene.sets, function(x) str_to_title(x))

gene.sets <- list()

# limma_rand2021.dt <- fread("/Users/rargelaguet/data/Rand2021_oskm_combinations_microarray/results/limma/1F2F/limma_output_all.csv.gz") %>%
#   .[feature%in%rownames(sce)]
# for (i in unique(limma_rand2021.dt$contrast)) {
#   gene.sets[[paste0("rand2021_",i)]] <- limma_rand2021.dt[contrast==i & pval_adj<0.01 & abs(logFC)>=1,feature]
# }

limma_rand2021.dt <- fread("/Users/rargelaguet/data/Rand2021_oskm_combinations_microarray/results/limma/pairwise_all/limma_results_vs_mock.txt.gz") %>%
  .[,feature:=str_to_title(feature)] %>%
  .[feature%in%rownames(sce)] %>% .[groupB!="iPSC"]

for (i in unique(limma_rand2021.dt$groupB)) {
  gene.sets[[paste0("rand2021_",i)]] <- limma_rand2021.dt[groupB==i & pvalue<0.05 & logFC>=1,feature] %>% head(n=1000)
}

sapply(gene.sets,length)

###############
## Run UCell ##
###############

sce <- ScoreSignatures_UCell(sce, features=gene.sets, maxRank=2500)

# Save
saveRDS(altExp(sce,'UCell'), file.path(io$outdir,"UCell_SummarizedExperiment_GRN.rds"))
