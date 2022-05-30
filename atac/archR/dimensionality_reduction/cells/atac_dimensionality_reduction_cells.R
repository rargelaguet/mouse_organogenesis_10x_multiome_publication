here::i_am("atac/archR/dimensionality_reduction/cells/atac_dimensionality_reduction_cells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(uwot))
suppressPackageStartupMessages(library(batchelor))

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
p$add_argument('--binarise', action="store_true",                                 help='Binarise ATAC matrix?')
p$add_argument('--scale_dims', action="store_true",                                 help='Scale latent dimensions?')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ndims',           type="integer",    default=30,                  help='Number of LSI dimensions')
p$add_argument('--batch_variable',  type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--batch_method',    type="character",  default="MNN",               help='Batch correctin method ("Harmony" or "MNN")')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$atac_matrix_file <- file.path(io$basedir,"processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
# args$atac_feature_stats <- file.path(io$basedir,"results/atac/archR/feature_stats/PeakMatrix/PeakMatrix_celltype_stats.txt.gz")
# args$stages <- "all"
# args$samples <- "all"
# args$matrix <- "PeakMatrix"
# args$nfeatures <- 25000
# args$remove_ExE_cells <- "False"
# args$batch_variable <- "sample"
# args$batch_method <- "MNN"
# args$binarise <- FALSE
# args$ndims <- 50
# args$scale_dims <- TRUE
# args$seed <- 42
# args$n_neighbors <- 25
# args$min_dist <- 0.50
# args$colour_by <- c("stage","celltype","nFrags_atac","pass_rnaQC")
# args$outdir <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/test")
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

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE & stage%in%args$stages & sample%in%args$samples] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

stopifnot(args$colour_by %in% colnames(sample_metadata))

# Sanity checks
if (length(args$batch_variable)>0) {
  stopifnot(args$batch_variable %in% colnames(sample_metadata))
}

# remove ExE cells
if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

# sanity checks
if (length(args$batch_variable)>0) {
  stopifnot(args$batch_variable%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch_variable]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch_variable))
    args$batch_variable <- NULL
    args$batch_method <- NULL
  } else {
    library(batchelor)
  }
}

# print stats
table(sample_metadata$stage)
table(sample_metadata$celltype)

#######################
## Fetch ATAC matrix ##
#######################

print(sprintf("Fetching single-cell ATAC %s...",args$matrix))

if (file.exists(args$atac_matrix_file)) {
  atac.se <- readRDS(args$atac_matrix_file)[,sample_metadata$cell]
} else {
  print(sprintf("%s does not exist. Loading matrix from the ArchR project...",args$atac_matrix_file))
  source(here::here("atac/archR/load_archR_project.R"))
  stopifnot(args$matrix %in% getAvailableMatrices(ArchRProject))
  ArchRProject <- ArchRProject[sample_metadata$cell]
  atac.se <- getMatrixFromProject(ArchRProject, useMatrix=args$matrix, binarize = FALSE)[,sample_metadata$cell]
  dim(atac.se)

  # Define feature names
  if (grepl("peak",tolower(args$matrix),ignore.case=T)) {
    row.ranges.dt <- rowRanges(atac.se) %>% as.data.table %>% 
      setnames("seqnames","chr") %>%
      .[,c("chr","start","end")] %>%
      .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
    rownames(atac.se) <- row.ranges.dt$idx
  } else if (grepl("gene",tolower(args$matrix),ignore.case=T)) {
    rownames(atac.se) <- rowData(atac.se)$name
  } else {
    stop("Matrix not recognised")
  }
}

#######################
## Feature selection ##
#######################

# Load feature stats
atac_featureStats.dt <- fread(args$atac_feature_stats)

# Define highly variable features using the pseudobulk estimates
atac_features <- atac_featureStats.dt %>% 
  setorder(-var_pseudobulk) %>% 
  head(n=args$nfeatures) %>% .$feature

stopifnot(atac_features%in%rownames(atac.se))

# if (args$test) {
#   print("Test mode activated, subsetting number of ATAC features...")
#   atac_features <- atac_features %>% head(n=500)
# }

#############################
## ATAC TFIDF normalisation ##
#############################

atac_tfidf.mtx <- tfidf(assay(atac.se[atac_features,]), method=1, scale.factor=1e4)

###########################
## Latent Semantic Index ##
###########################

svd <- irlba::irlba(atac_tfidf.mtx, args$ndims, args$ndims)
svdDiag <- matrix(0, nrow=args$ndims, ncol=args$ndims)
diag(svdDiag) <- svd$d
lsi.mtx <- t(svdDiag %*% t(svd$v))
rownames(lsi.mtx) <- colnames(atac.se)
colnames(lsi.mtx) <- paste0("LSI",seq_len(ncol(lsi.mtx)))

rm(svd,svdDiag); gc(reset=T)

###############
## Parse LSI ##
###############

# (Optional) Remove dimensions that correlate with nFrags
if (opts$remove_dim_cor_seq_depth) {
  lsi.mtx <- lsi.mtx[,which(abs(cor(lsi.mtx,sample_metadata$nFrags_atac))<=0.75)]
}

# (Optional) scale dimensions to Z scores and remove outliers
# lsi.mtx <- sweep(lsi.mtx - rowMeans(lsi.mtx), 1, matrixStats::rowSds(lsi.mtx),`/`)
if (args$scale_dims) {
  lsi.mtx <- sweep(lsi.mtx - colMeans(lsi.mtx), 2, matrixStats::colSds(lsi.mtx),`/`)
  lsi.mtx[lsi.mtx>2] <- 2
  lsi.mtx[lsi.mtx<(-2)] <- (-2)
}

######################
## Batch correction ##
######################

if (length(args$batch_variable)>0) {
  print(sprintf("Applying %s batch correction for variable: %s", args$batch_method, args$batch_variable))
  # Harmony
  if (args$batch_method=="Harmony") {
    stop("Not implemented")
    # (...)
    # library(harmony)
    # harmonyParams <- list(...)
    # harmonyParams$data_mat <- getReducedDims(
    #   ArchRProj = ArchRProj, 
    #   reducedDims = reducedDims, 
    #   dimsToUse = dimsToUse, 
    #   scaleDims = scaleDims, 
    #   corCutOff = corCutOff
    # )
    # harmonyParams$verbose <- verbose
    # harmonyParams$meta_data <- data.frame(getCellColData(
    #   ArchRProj = ArchRProj, 
    #   select = groupBy)[rownames(harmonyParams$data_mat), , drop = FALSE])
    # harmonyParams$do_pca <- FALSE
    # harmonyParams$vars_use <- groupBy
    # harmonyParams$plot_convergence <- FALSE
    
    # lsi.dt <- getReducedDims(ArchRProject, "IterativeLSI_Harmony") %>% round(3) %>% 
    #   as.data.table(keep.rownames = T) %>% setnames("rn","cell")

  } else if (args$batch_method=="MNN") {
      
    stopifnot(args$batch_variable=="sample")
    
    # Define stage and sample order
    stopifnot(sample_metadata$cell==rownames(lsi.mtx))
    timepoints <- sample_metadata$stage
    timepoint_order <- opts$stages[opts$stages%in%timepoints]
    samples <- sample_metadata$sample
    sample_order <- opts$samples[opts$samples%in%samples]
    
    lsi_list    <- lapply(unique(timepoints), function(i){
      sub_pc   <- lsi.mtx[timepoints == i, , drop = FALSE]
      sub_samp <- samples[timepoints == i]
      list     <- lapply(unique(sub_samp), function(j){ sub_pc[sub_samp == j, , drop = FALSE]})
      names(list) <- unique(sub_samp)
      return(list)
    })
    names(lsi_list) <- unique(timepoints)
    
    #arrange to match timepoint order
    lsi_list <- lsi_list[order(match(names(lsi_list), timepoint_order))]
    lsi_list <- lapply(lsi_list, function(x){ x[order(match(names(x), sample_order))]})
    
    #perform corrections within stages
    correct_list <- lapply(lsi_list, function(x){
      if(length(x) > 1){
        return(do.call(reducedMNN, x)$corrected)
      } else {
        return(x[[1]])
      }
    })
    
    # perform correction over stages
    lsi.mtx <- reducedMNN(correct_list, merge.order=1:length(correct_list))$corrected 
    rm(correct_list,lsi_list)
    
  } else {
    stop("Batch correction method not recognised")
  }
}

# Save LSI coordinates
lsi.dt <- lsi.mtx %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(lsi.dt, sprintf("%s/lsi_nfeatures%d_ndims%d.txt.gz",args$outdir, args$nfeatures, args$ndims))

##########
## UMAP ##
##########

pt.size <- ifelse(nrow(lsi.mtx)>=1e4,0.8,1.2)

# i <- args$n_neighbors[1]; j <- args$min_dist[1]
for (i in args$n_neighbors) {
  for (j in args$min_dist) {
    
    # Define the latent space to run UMAP on
    if (length(opts$batch.correction)>0) {
      if (args$batch_method=="Harmony") {
        dimred <- "IterativeLSI_Harmony"
      } else if  (args$batch_method=="MNN") {
        dimred <- "IterativeLSI_MNN"
      } else {
        stop("Batch correction method not recognised")
      }
    } else {
      dimred <- "IterativeLSI"
    }
    
    # Run UMAP
    set.seed(args$seed)
    umap_embedding.mtx <- umap(lsi.mtx, n_neighbors=i, min_dist=j, metric="cosine", fast_sgd = TRUE) %>% round(2)
    rownames(umap_embedding.mtx) <- rownames(lsi.mtx)
    
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
      merge(sample_metadata,by="cell")
    
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

# to.plot <- umap.dt %>% merge(sample_metadata,by="cell")

# # to.plot[,foo:=F]
# # to.plot[pass_rnaQC==F & nFrags_atac_log10<=4,foo:=TRUE]
# # to.plot[,foo:=nFrags_atac<=10000]

# ggplot(to.plot, aes_string(x="umap1", y="umap2", fill="PromoterRatio_atac")) +
#   geom_point(size=1.5, shape=21, stroke=0.05) +
#   scale_fill_gradientn(colours = terrain.colors(10)) +
#   theme_classic() +
#   ggplot_theme_NoAxes()
