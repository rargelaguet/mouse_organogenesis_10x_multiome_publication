here::i_am("rna/dimensionality_reduction/dimensionality_reduction_sce.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(Seurat))

# TO-DO: BATCH CORRECTION, EXTRACT DATA FROM DIMRED AND GGPLOT2

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--seurat',             type="character",                               help='Seurat object file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--features',        type="integer",    default=1000,                help='Number of features')
p$add_argument('--npcs',            type="integer",    default=30,                  help='Number of PCs')
p$add_argument('--vars_to_regress', type="character",                nargs='+',     help='Metadata columns to regress out')
p$add_argument('--n_neighbors',     type="integer",    default=30,     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  default="celltype.mapped",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output file')
p$add_argument('--remove_ExE_cells',       type="character",  default="False",  help='Remove ExE cells? ("True"/"False")')
# p$add_argument('--batch_correction', type="character",                               help='Metadata column to apply batch correction on')
# p$add_argument('--SCTransform', action="store_true",                                 help='Remove ExE cells?')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args$seurat <- io$rna.seurat
# args$metadata <- "/bi/group/reik/ricard/data/gastrulation_histones/results/rna/mapping/sample_metadata_after_mapping.txt.gz"
# args$stages <- "E8.75"
# args$features <- 1000
# args$npcs <- 30
# args$colour_by <- c("celltype.mapped","celltype.mapped_seurat","sample")
# args$vars_to_regress <- c("nFeature_RNA","mitochondrial_percent_RNA")
# args$n_neighbors <- 25
# args$min_dist <- 0.5
# args$seed <- 42
# args$outdir <- "/bi/group/reik/ricard/data/gastrulation_histones/results/rna/dimensionality_reduction/seurat"
# args$remove_ExE_cells <- TRUE
## END TEST ##

# if (isTRUE(args$test)) print("Test mode activated...")

#####################
## Parse arguments ##
#####################

if (args$stages[1]=="all") {
  args$stages <- opts$stages
} else {
  stopifnot(args$stages%in%opts$stages)
}

if (args$remove_ExE_cells=="True") {
  args$remove_ExE_cells <- TRUE
} else if (args$remove_ExE_cells=="False") {
  args$remove_ExE_cells <- FALSE 
} else {
  stop('remove_ExE_cells should be "True" or "False"')
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  # .[nFeature_RNA>=2000] %>% .[stage=="E8.7",stage:="E8.75"] %>%    # temporary
  # .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples]
  .[pass_rnaQC==TRUE & doublet_call==FALSE & stage%in%args$stages]

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

table(sample_metadata$stage)
table(sample_metadata$celltype.mapped)

###################
## Sanity checks ##
###################

stopifnot(args$colour_by %in% colnames(sample_metadata))
# stopifnot(unique(sample_metadata$celltype.mapped) %in% names(opts$celltype.colors))

# if (length(args$batch.correction)>0) {
#   stopifnot(args$batch.correction%in%colnames(sample_metadata))
#   if (length(unique(sample_metadata[[args$batch.correction]]))==1) {
#     message(sprintf("There is a single level for %s, no batch correction applied",args$batch.correction))
#     args$batch.correction <- NULL
#   } else {
#     library(batchelor)
#   }
# }

if (length(args$vars_to_regress)>0) {
  stopifnot(args$vars_to_regress%in%colnames(sample_metadata))
}


###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
seurat <- load_Seurat(
  file = args$seurat, 
  cells = sample_metadata$cell,
  normalise = TRUE, scale = FALSE, 
  remove_non_expressed_genes = TRUE
)
dim(seurat)

# Update sample metadata
foo <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>% .[colnames(seurat),]
stopifnot(colnames(seurat) == rownames(foo))
seurat@meta.data <- foo

#################
## SCTransform ##
#################

#   variable.features.n = args$features, 
#   # vars_to_regress = c("nFeature_RNA","mitochondrial_percent_RNA"),
#   vars_to_regress = args$vars_to_regress,
#   verbose = FALSE
# )

#######################
## Feature selection ##
#######################

seurat <- FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = args$features)
# seurat <- FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = 1000, assay = "SCT")
# head(seurat@assays$RNA@var.features)

###########################################
## Scale data and regress out covariates ##
###########################################

# seurat <- ScaleData(seurat, features=var.genes, vars_to_regress=c("nFeature_RNA","mitochondrial_percent_RNA"))
seurat <- ScaleData(seurat, 
  features = VariableFeatures(seurat), 
  vars.to.regress = args$vars_to_regress, 
  verbose = FALSE
)

#########
## PCA ##
#########

seurat <- RunPCA(seurat, features = VariableFeatures(seurat), npcs = args$npcs, verbose = FALSE)

# Save PCA coordinates
pca.dt <- seurat@reductions[["pca"]]@cell.embeddings %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(pca.dt, sprintf("%s/pca_features%d_pcs%d.txt.gz",args$outdir, args$features, args$npcs))


##########
## UMAP ##
##########

# Run
# set.seed(args$seed)
seurat <- RunUMAP(seurat, 
  dims = 1:args$npcs,
  reduction = "pca",
  n.neighbors = args$n_neighbors,
  min.dist = args$min_dist
)

# Fetch UMAP coordinates
umap.dt <- seurat@reductions[["umap"]]@cell.embeddings %>% round(3) %>% as.data.table %>% 
  .[,cell:=colnames(seurat)] %>%
  setnames(c("UMAP1","UMAP2","cell")) %>% 
  .[,c("cell","UMAP1","UMAP2")]

# Save UMAP coordinates
fwrite(umap.dt, sprintf("%s/umap_features%d_pcs%d_neigh%d_dist%s.txt.gz",args$outdir, args$features, args$npcs, args$n_neighbors, args$min_dist))

##########
## Plot ##
##########

pt.size <- ifelse(ncol(seurat)>=1e4,0.5,0.75)

for (i in args$colour_by) {

  Idents(seurat) <- i

  p <- DimPlot(seurat, label = FALSE, reduction = 'umap', pt.size = pt.size) + 
    NoAxes()

  # Define colormap
  if (is.numeric(seurat[[i]])) {
    p <- p + scale_fill_gradientn(colours = terrain.colors(10))
  }
  if (grepl("celltype",i)) {
    p <- p + scale_colour_manual(values=opts$celltype.colors) +
      theme(
        legend.position="none",
        legend.title=element_blank()
      )
  }
  if (grepl("stage",i)) {
    p <- p + scale_colour_manual(values=opts$stage.colors) +
      guides(fill = guide_legend(override.aes = list(size=3))) +
      theme(
        legend.position="right",
        legend.title=element_blank()
      )
  }

  # Save UMAP plot
  outfile <- file.path(args$outdir,sprintf("umap_features%d_pcs%d_neigh%d_dist%s_%s.pdf", args$features, args$npcs, args$n_neighbors, args$min_dist, i))
  pdf(outfile, width=7, height=5)
  print(p)
  dev.off()
}

