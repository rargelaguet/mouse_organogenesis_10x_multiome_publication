here::i_am("rna/dimensionality_reduction/metacells/dimensionality_reduction_sce_metacells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--features',        type="integer",    default=1000,                help='Number of features')
p$add_argument('--npcs',            type="integer",    default=30,                  help='Number of PCs')
p$add_argument('--n_neighbors',     type="integer",    default=30,     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  default="celltype",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--remove_ExE_cells',       type="character",  default="False",  help='Remove ExE cells? ("True"/"False")')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output file')
# p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--vars_to_regress', type="character",                nargs='+',     help='Metadata columns to regress out')
p$add_argument('--batch_correction',type="character",                               help='Metadata column to apply batch correction on')
# p$add_argument('--test',      action = "store_true",                       help='Testing mode')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$sce <- file.path(io$basedir,"results/rna/metacells/all_cells/SingleCellExperiment_metacells.rds") # io$rna.sce
# args$metadata <- file.path(io$basedir,"results/rna/metacells/all_cells/metacells_metadata.txt.gz")
# args$stages <- c("E8.5","E8.75") # "all" # "E7.75"
# args$features <- 2500
# args$npcs <- 50
# args$colour_by <- c("celltype")
# args$vars_to_regress <- NULL # c("nFeature_RNA","mitochondrial_percent_RNA")
# args$batch_correction <- "sample"
# args$remove_ExE_cells <- "False"
# args$n_neighbors <- 25
# args$min_dist <- 0.5
# args$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction/metacells")
## END TEST ##

#####################
## Define settings ##
#####################

# if (isTRUE(args$test)) print("Test mode activated...")
dir.create(args$outdir, showWarnings = F)

#####################
## Parse arguments ##
#####################

# Options
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

sample_metadata <- fread(args$metadata) %>% .[genotype=="WT"]

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

table(sample_metadata$stage)
table(sample_metadata$celltype)
table(sample_metadata$sample)

###################
## Sanity checks ##
###################

stopifnot(args$colour_by %in% colnames(sample_metadata))
# stopifnot(unique(sample_metadata$celltype) %in% names(opts$celltype.colors))

if (length(args$batch_correction)>0) {
  stopifnot(args$batch_correction%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch_correction]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch_correction))
    args$batch_correction <- NULL
  } else {
    library(batchelor)
  }
}

if (length(args$vars_to_regress)>0) {
  stopifnot(args$vars_to_regress%in%colnames(sample_metadata))
}


###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- readRDS(args$sce)[,sample_metadata$metacell]

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("metacell") %>% DataFrame

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
# if (length(args$batch_correction)>0) {
#   decomp <- modelGeneVar(sce, block=colData(sce)[[args$batch_correction]])
# } else {
#   decomp <- modelGeneVar(sce)
# }
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=args$features) %>% rownames

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]

############################
## Regress out covariates ##
############################

if (length(args$vars_to_regress)>0) {
  print(sprintf("Regressing out variables: %s", paste(args$vars_to_regress,collapse=" ")))
  logcounts(sce_filt) <- RegressOutMatrix(
    mtx = logcounts(sce_filt),
    covariates = colData(sce_filt)[,args$vars_to_regress,drop=F]
  )
}

############################
## PCA + Batch correction ##
############################

# outfile <- sprintf("%s/%s_pca_features%d_pcs%d.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$features, args$npcs)
if (length(args$batch_correction)>0) {
  suppressPackageStartupMessages(library(batchelor))
  print(sprintf("Applying MNN batch correction for variable: %s", args$batch_correction))
  pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch_correction]], d = args$npcs)
  pca.corrected <- reducedMNN(pca)$corrected
  colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
  reducedDim(sce_filt, "PCA") <- pca.corrected[colnames(sce),]
} else {
  sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)
}

# Save PCA coordinates
pca.dt <- reducedDim(sce_filt,"PCA") %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(pca.dt, sprintf("%s/pca_features%d_pcs%d.txt.gz",args$outdir, args$features, args$npcs))

##########
## UMAP ##
##########

# Run
set.seed(args$seed)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = args$n_neighbors, min_dist = args$min_dist)

# Fetch UMAP coordinates
umap.dt <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  setnames(c("UMAP1","UMAP2","cell"))

# Save UMAP coordinates
fwrite(umap.dt, sprintf("%s/umap_features%d_pcs%d_neigh%d_dist%s.txt.gz",args$outdir, args$features, args$npcs, args$n_neighbors, args$min_dist))

##########
## Plot ##
##########

pt.size <- ifelse(ncol(sce)>=1e4,1.2,2)

for (i in args$colour_by) {

  to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
    .[,metacell:=colnames(sce_filt)] %>%
    merge(sample_metadata, by="metacell")

  if (is.numeric(to.plot[[i]])) {
    if (max(to.plot[[i]],na.rm=T) - min(to.plot[[i]],na.rm=T) > 1000) {
      to.plot[[i]] <- log10(to.plot[[i]]+1)
      to.plot %>% setnames(i,paste0(i,"_log10")); i <- paste0(i,"_log10")
    }
  }
  
  p <- ggplot(to.plot, aes_string(x="V1", y="V2", fill=i)) +
    geom_point(size=2, shape=21, stroke=0.05) +
    theme_classic() +
    ggplot_theme_NoAxes()
  
  # Define colormap
  if (is.numeric(to.plot[[i]])) {
    p <- p + scale_fill_gradientn(colours = terrain.colors(10))
  }
  if (grepl("celltype",i)) {
    p <- p + scale_fill_manual(values=opts$celltype.colors) +
      theme(
        legend.position="none",
        legend.title=element_blank()
      )
  }
  if (grepl("stage",i)) {
    p <- p + scale_fill_manual(values=opts$stage.colors) +
      guides(fill = guide_legend(override.aes = list(size=3))) +
      theme(
        legend.position="right",
        legend.title=element_blank()
      )
  }

  # Save UMAP plot
  outfile <- file.path(args$outdir,sprintf("umap_features%d_pcs%d_neigh%d_dist%s_%s.pdf",args$features, args$npcs, args$n_neighbors, args$min_dist, i))
  pdf(outfile, width=7, height=5)
  print(p)
  dev.off()
}



##########
## Test ##
##########

# to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
#   .[,metacell:=colnames(sce_filt)]# %>%
#   
# to.plot$expr <- logcounts(sce)["Pax7",]*logcounts(sce)["T",]
# to.plot[expr>=30,expr:=30]
# 
# to.plot <- to.plot %>% merge(sample_metadata, by="metacell")
# 
# ggplot(to.plot, aes_string(x="V1", y="V2", fill="expr")) +
#   geom_point(size=2, shape=21, stroke=0.05) +
#   # viridis::scale_fill_viridis() +
#   scale_fill_gradient(low = "gray80", high = "purple") +
#   theme_classic() +
#   ggplot_theme_NoAxes()
# 
# 
# ggplot(to.plot, aes_string(x="celltype", y="expr", fill="celltype")) +
#   geom_boxplot(outlier.shape=NA) +
#   labs(x="", y="T/Pax7 coexpression") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   theme_classic() +
#   guides(x = guide_axis(angle = 90)) +
#   theme(
#     legend.position = "none",
#     axis.text.y = element_text(color="black", size=rel(1.0)),
#     axis.text.x = element_text(color="black", size=rel(0.8)),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )
