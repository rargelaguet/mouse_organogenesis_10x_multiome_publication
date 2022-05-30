here::i_am("atac/archR/dimensionality_reduction/metacells/atac_dimensionality_reduction_metacells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(uwot))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--matrix',          type="character",   help='Matrix to use')
p$add_argument('--atac_matrix_file',          type="character",  help='Matrix file')
p$add_argument('--atac_feature_stats',          type="character",   help='Feature stats')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--samples',       type="character",  default="all",  nargs='+',  help='Samples to plot')
p$add_argument('--remove_ExE_cells',       type="character",  default="False",  help='Remove ExE cells? ("True"/"False")')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ndims',           type="integer",    default=30,                  help='Number of LSI dimensions')
p$add_argument('--batch_variable',  type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--batch_method',    type="character",  default="MNN",               help='Batch correctin method ("Harmony" or "MNN")')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$metadata <- file.path(io$basedir,"results/atac/archR/metacells/all_cells/PeakMatrix/metacells_metadata.txt.gz")
# args$matrix <- "PeakMatrix"
# args$atac_feature_stats <- file.path(io$basedir,sprintf("results/atac/archR/feature_stats/%s/%s_celltype_stats.txt.gz",args$matrix,args$matrix))
# args$atac_matrix_file <- file.path(io$basedir,sprintf("results/atac/archR/metacells/all_cells/%s/%s_summarized_experiment_metacells.rds",args$matrix,args$matrix))
# args$stages <- "E8.75" # "all"
# args$samples <- "all"
# args$remove_ExE_cells <- "True"
# args$nfeatures <- 15000
# args$batch_variable <- "None"
# args$batch_method <- "MNN"
# args$binarise <- FALSE
# args$ndims <- 50
# args$n_neighbors <- 25
# args$min_dist <- 0.50
# args$colour_by <- c("celltype")
# args$outdir <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/metacells")
## END TEST ##

#####################
## Parse arguments ##
#####################

if (args$stages[1]=="all") {
  args$stages <- opts$stages
} else {
  stopifnot(args$stages%in%opts$stages)
}

if (args$samples[1]=="all") {
  args$samples <- opts$samples
} else {
  stopifnot(args$samples%in%opts$samples)
}


if (args$remove_ExE_cells=="True") {
  args$remove_ExE_cells <- TRUE
} else if (args$remove_ExE_cells=="False") {
  args$remove_ExE_cells <- FALSE 
} else {
  stop('remove_ExE_cells should be "True" or "False"')
}

if (args$batch_variable=="None") {
  args$batch_variable <- NULL
  args$batch_method <- NULL
}

print(args)

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

# Options
opts$remove_dim_cor_seq_depth <- TRUE

##########################
## Load sample metadata ##
##########################

metacell_metadata.dt <- fread(args$metadata) %>%
  .[stage%in%args$stages & genotype=="WT"] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

stopifnot(args$colour_by %in% colnames(metacell_metadata.dt))

# Sanity checks
if (length(args$batch_variable)>0) {
  stopifnot(args$batch_variable %in% colnames(metacell_metadata.dt))
}

# remove ExE cells
if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  metacell_metadata.dt <- metacell_metadata.dt %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

# sanity checks
if (length(args$batch_variable)>0) {
  stopifnot(args$batch_variable%in%colnames(metacell_metadata.dt))
  if (length(unique(metacell_metadata.dt[[args$batch_variable]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch_variable))
    args$batch_variable <- NULL
    args$batch_method <- NULL
  } else {
    library(batchelor)
  }
}
# print stats
table(metacell_metadata.dt$stage)
table(metacell_metadata.dt$celltype)

#######################
## Fetch ATAC matrix ##
#######################

print(sprintf("Fetching metacell ATAC %s...",args$matrix))

atac.se <- readRDS(args$atac_matrix_file)[,metacell_metadata.dt$metacell]

# Normalise
assay(atac.se) <- log(1e6*(sweep(assay(atac.se),2,colSums(assay(atac.se),na.rm=T),"/"))+1)

#######################
## Feature selection ##
#######################

# Load feature stats
atac_featureStats.dt <- fread(args$atac_feature_stats)

# Define highly variable features using the pseudobulk estimates
atac_features <- atac_featureStats.dt %>% 
  setorder(-var_metacells) %>% 
  head(n=args$nfeatures) %>% .$feature

stopifnot(atac_features%in%rownames(atac.se))

#########
## PCA ##
#########

pca.mtx <- irlba::prcomp_irlba(t(assay(atac.se)), n=args$ndims)$x
rownames(pca.mtx) <- colnames(atac.se)
colnames(pca.mtx) <- paste0("PC",seq_len(ncol(pca.mtx)))

###############
## Parse LSI ##
###############

# (Optional) Remove dimensions that correlate with nFrags
if (opts$remove_dim_cor_seq_depth) {
  pca.mtx <- pca.mtx[,which(cor(pca.mtx,metacell_metadata.dt$nFrags_atac)<=0.75)]
}

# (Optional) scale dimensions, Z scores
# pca.mtx <- sweep(pca.mtx - rowMeans(pca.mtx), 1, matrixStats::rowSds(pca.mtx),`/`)
# if (args$scale_dims) {
#   pca.mtx <- sweep(pca.mtx - colMeans(pca.mtx), 2, matrixStats::colSds(pca.mtx),`/`)
#   pca.mtx[pca.mtx>2] <- 2
#   pca.mtx[pca.mtx<(-2)] <- (-2)
# }

######################
## Batch correction ##
######################

if (length(args$batch_variable)>0) {
  print(sprintf("Applying %s batch correction for variable: %s", args$batch_method, args$batch_variable))
  if (args$batch_method=="Harmony") {
    stop("Not implemented")
  } else if (args$batch_method=="MNN") {
    library(batchelor)
    stopifnot(rownames(pca.mtx)==metacell_metadata.dt$cell)
    pca.mtx <- reducedMNN(pca.mtx, batch=metacell_metadata.dt[[args$batch_variable]])$corrected
    colnames(pca.mtx) <- paste0("LSI",seq_len(ncol(pca.mtx)))
    pca.mtx <- pca.mtx[colnames(atac.se),]
  } else {
    stop("Batch correction method not recognised")
  }
}


# Save LSI coordinates
# outfile <- sprintf("%s/lsi_%s_nfeatures%d_ndims%d.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims)
 # outfile <- sprintf("%s/lsi_nfeatures%d_dims%d_%sbatchcorrection_by_%s.txt.gz",args$outdir, args$nfeatures, args$ndims, args$batch_method, paste(args$batch_variable,collapse="-"))
pca.dt <- pca.mtx %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
outfile <- sprintf("%s/pca_nfeatures%d_ndims%d.txt.gz",args$outdir, args$nfeatures, args$ndims)
fwrite(pca.dt, outfile)

##########
## UMAP ##
##########

pt.size <- ifelse(nrow(pca.mtx)>=1e3,1.5,2.5)

# i <- args$n_neighbors[1]; j <- args$min_dist[1]
for (i in args$n_neighbors) {
  for (j in args$min_dist) {

    # Run UMAP
    set.seed(42)
    umap_embedding.mtx <- umap(pca.mtx, n_neighbors=i, min_dist=j, metric="cosine", fast_sgd = TRUE) %>% round(2)
    rownames(umap_embedding.mtx) <- rownames(pca.mtx)
    
    # Fetch UMAP coordinates
    umap.dt <- umap_embedding.mtx %>%
      as.data.table(keep.rownames = T) %>%
      setnames(c("cell","umap1","umap2"))
    
    # Save UMAP coordinates
    # outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_neigh%d_dist%s.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims, i, j)
    outfile <- sprintf("%s/umap_nfeatures%d_ndims%d.txt.gz",args$outdir, args$nfeatures, args$ndims)
    fwrite(umap.dt, outfile)

    # Plot
    to.plot <- umap.dt %>%
      merge(metacell_metadata.dt,by="cell")
    
    # k <- "celltype"
    for (k in args$colour_by) {

      # log10 large numeric values
      if (is.numeric(to.plot[[k]])) {
        if (max(to.plot[[k]],na.rm=T) - min(to.plot[[k]],na.rm=T) > 1000) {
          to.plot[[k]] <- log10(to.plot[[k]]+1)
          to.plot %>% setnames(k,paste0(k,"_log10")); k <- paste0(k,"_log10")
        }
      }
      
      p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", fill=k)) +
        geom_point(size=pt.size, shape=21, stroke=0.05) +
        # ggrastr::geom_point_rast(size=1.5, shape=21, stroke=0.05) +  # DOES NOT WORK IN THE CLUSTER
        theme_classic() +
        ggplot_theme_NoAxes()
      
      # Define colormap
      if (is.numeric(to.plot[[j]])) {
        p <- p + scale_fill_gradientn(colours = terrain.colors(10))
      }

      if (grepl("celltype",k)) {
        p <- p + scale_fill_manual(values=opts$celltype.colors) +
          theme(
            legend.position="none",
            legend.title=element_blank()
          )
      }
      if (grepl("stage",i)) {
        p <- p + scale_fill_manual(values=opts$stage.colors) +
          theme(
            legend.position="none",
            legend.title=element_blank()
          )
      }
      
      # Save UMAP plot
      # outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_neigh%d_dist%s_%s.pdf",args$outdir, args$matrix, args$nfeatures, args$ndims, i, j, k)
      outfile <- sprintf("%s/umap_nfeatures%d_ndims%d_%s.pdf",args$outdir, args$nfeatures, args$ndims, k)
      pdf(outfile, width=7, height=5)
      print(p)
      dev.off()
    }
    
  }
}


##########
## TEST ##
##########

# to.plot <- umap.dt %>% merge(metacell_metadata.dt,by="cell")

# # to.plot[,foo:=F]
# # to.plot[pass_rnaQC==F & nFrags_atac_log10<=4,foo:=TRUE]
# # to.plot[,foo:=nFrags_atac<=10000]

# ggplot(to.plot, aes_string(x="umap1", y="umap2", fill="PromoterRatio_atac")) +
#   geom_point(size=1.5, shape=21, stroke=0.05) +
#   scale_fill_gradientn(colours = terrain.colors(10)) +
#   theme_classic() +
#   ggplot_theme_NoAxes()
