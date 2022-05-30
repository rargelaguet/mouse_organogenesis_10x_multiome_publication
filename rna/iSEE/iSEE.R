source(here::here("settings.R"))
source(here::here("utils.R"))

library(scran)
library(scater)
library(shiny)
library(iSEE)

#####################
## Define settings ##
#####################

opts$samples <- c(
  # "E7.5_rep1",
  # "E7.5_rep2",
  "E7.75_rep1"
  # "E8.0_rep1",
  # "E8.0_rep2",
  # "E8.5_rep1",
  # "E8.5_rep2",
  # "E8.75_rep1",
  # "E8.75_rep2",
  # "E8.5_CRISPR_T_KO",
  # "E8.5_CRISPR_T_WT"
)

opts$remove_ExE_cells <- FALSE

opts$vars_to_regress <- NULL # c("nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA")

##########################
## Load sample metadata ##
##########################

# io$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz") # io$metadata
io$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
io$rna.sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")

sample_metadata <- fread(io$metadata) %>%
  # .[ribosomal_percent_RNA<=5] %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples]

if (opts$remove_ExE_cells) {
  sample_metadata <- sample_metadata %>% .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

table(sample_metadata$sample)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

# Filter out some genes manually
sce <- sce[!grepl("Rpl|Rps|mt-",rownames(sce)),]

# Select highly variable genes
decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=2500) %>% rownames
sce_filt <- sce[hvgs,]

############################
## Regress out covariates ##
############################

if (length(opts$vars_to_regress)>0) {
  print(sprintf("Regressing out variables: %s", paste(opts$vars_to_regress,collapse=" ")))
  logcounts(sce_filt) <- RegressOutMatrix(
    mtx = logcounts(sce_filt),
    covariates = colData(sce_filt)[,opts$vars_to_regress,drop=F]
  )
}

##############################
## Dimensionality reduction ##
##############################

# PCA
# npcs <- 30
# data <- scale(t(logcounts(sce_filt)), center = T, scale = F)
# reducedDim(sce_filt, "PCA") <- irlba::prcomp_irlba(data, n=npcs)$x#[,1:npcs]
sce_filt <- runPCA(sce_filt, ncomponents = 30, ntop = 2500)  
reducedDim(sce,"PCA") <- reducedDim(sce_filt,"PCA")

# UMAP
set.seed(42)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 30, min_dist = 0.3)
reducedDim(sce,"UMAP") <- reducedDim(sce_filt,"UMAP")

plotUMAP(sce_filt, colour_by="celltype.mapped") + 
  scale_color_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none"
  )

plotUMAP(sce_filt, colour_by="ribosomal_percent_RNA")
plotUMAP(sce_filt, colour_by="mitochondrial_percent_RNA")
plotUMAP(sce_filt, colour_by="nFeature_RNA")

sce_filt$nFrags_atac <- log10(sce_filt$nFrags_atac)
plotUMAP(sce_filt, colour_by="nFrags_atac")

plot(sce_filt$nFeature_RNA, sce_filt$mitochondrial_percent_RNA)

#######################
## Define color maps ##
#######################

celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(sce$celltype.mapped)]
stopifnot(unique(sce$celltype.mapped) %in% names(celltype.colors))
stopifnot(names(celltype.colors)%in%unique(sce$celltype.mapped))
sce$celltype.mapped <- factor(sce$celltype.mapped, levels=names(celltype.colors))

celltype_color_fun <- function(n){
  return(celltype.colors)
}

categorical_color_fun <- function(n){
  return(RColorBrewer::brewer.pal(n, "Set2"))
}

# Define color maps
ecm <- ExperimentColorMap(
  # List of colormaps for assays.
  # assays = list(
  #   counts = viridis::viridis,
  #   cufflinks_fpkm = fpkm_color_fun
  # ),
  colData = list(
    celltype.mapped = celltype_color_fun
  ),
  # Colormaps applied to all undefined continuous assays
  all_continuous = list(
    assays = viridis::viridis
  ),
  # Colormaps applied to all undefined categorical assays
  all_discrete = list(
    assays = categorical_color_fun
  )
  # Colormap applied to all undefined categorical covariates.
  # global_discrete <- list()
  # Colormap applied to all undefined continuous covariates.
  # global_continuous <- list()
)

#########################
## Define iSEE options ##
#########################

# sce_filt <- registerAppOptions(sce_filt, color.maxlevels=40)
# getAppOption("color.maxlevels", sce_filt)

##############
## Run iSEE ##
##############

app <- iSEE(sce, colormap = ecm)
runApp(app)

##############################
## Load precomputed session ##
##############################

tmp <- readRDS("/Users/argelagr/Downloads/iSEE_memory.rds")
app <- iSEE(se = tmp$se, initial = tmp$memory, colormap = ecm)
runApp(app)
