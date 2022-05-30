here::i_am("atac/archR/dimensionality_reduction/cells/archR_dimensionality_reduction.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
# p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--remove_ExE_cells', action="store_true",                                 help='Remove ExE cells?')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ndims',           type="integer",    default=30,                  help='Number of LSI dimensions')
p$add_argument('--batch.variable',  type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--batch.method',    type="character",  default="MNN",               help='Batch correctin method ("Harmony" or "MNN")')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$metadata <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
# args$nfeatures <- 15000
# args$matrix <- "PeakMatrix"
# args$ndims <- 25
# args$seed <- 42
# args$n_neighbors <- 25
# args$min_dist <- 0.3
# args$colour_by <- c("sample","cluster")
# args$vars_to_regress <- c("nFeature_RNA","mitochondrial_percent_RNA")
# args$outdir <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results/atac/archR/dimensionality_reduction"
## END TEST ##


#####################
## Define settings ##
#####################

# Options
opts$lsi.iterations = 2
opts$lsi.cluster.resolution = 2

if (args$stages[1]=="all") {
  args$stages <- opts$stages
} else {
  stopifnot(args$stages%in%opts$stages)
}


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  # .[pass_atacQC==TRUE & doublet_call==FALSE & sample%in%args$samples] %>%
  # .[nFrags_atac>=3500 & TSSEnrichment_atac>=9,pass_atacQC:=TRUE] %>% .[stage=="E8.7",stage:="E8.75"] %>%    # temporary
  .[pass_atacQC==TRUE & doublet_call==FALSE & stage%in%args$stages] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

stopifnot(args$colour_by %in% colnames(sample_metadata))

if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

table(sample_metadata$stage)
table(sample_metadata$celltype)

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]
# ArchRProject.filt@sampleColData <- ArchRProject.filt@sampleColData[args$samples,,drop=F]

###################
## Sanity checks ##
###################

stopifnot(args$matrix %in% getAvailableMatrices(ArchRProject))

if (length(args$batch.variable)>0) {
  stopifnot(args$batch.variable%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch.variable]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch.variable))
    args$batch.variable <- NULL
  } else {
    library(batchelor)
  }
}

if (length(args$vars_to_regress)>0) {
  stopifnot(args$vars_to_regress%in%colnames(sample_metadata))
}

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(all(rownames(sample_metadata.to.archr) == rownames(getCellColData(ArchRProject.filt))))
for (i in args$colour_by) {
  ArchRProject.filt <- addCellColData(
    ArchRProject.filt,
    data = sample_metadata.to.archr[[i]],
    name = i,
    cells = rownames(sample_metadata.to.archr),
    force = TRUE
  )
  print(table(getCellColData(ArchRProject.filt,i)[[1]]))
}


#######################
## Feature selection ##
#######################

# if (args$matrix=="PeakMatrix_filt") {
  
#   # Load peak variability estimates
#   peak.variability.dt <- fread(io$archR.peak.variability)
  
#   # Define highly variable peaks
#   peaks <- peak.variability.dt %>% 
#     .[,peak:=sub("_",":",peak)] %>% .[,peak:=sub("_","-",peak)] %>%
#     setorder(-variance_pseudobulk) %>% 
#     head(n=args$nfeatures) %>% 
#     .$peak
  
#   # Subset peaks in the ArchR object
#   names(ArchRProject.filt@peakSet) <- sprintf("%s:%s-%s",seqnames(ArchRProject.filt@peakSet), start(ArchRProject.filt@peakSet), end(ArchRProject.filt@peakSet))
#   ArchRProject.filt@peakSet[peaks,] <- ArchRProject.filt@peakSet[peaks,]
  
#   ArchRProject.filt <- addFeatureMatrix(
#     input = ArchRProject.filt,
#     features = ArchRProject.filt@peakSet,
#     matrixName = "PeakMatrix_filt",
#     binarize = TRUE
#   )
  
# }

###########################
## Latent Semantic Index ##
###########################

# Iterative LSI: two iterations
ArchRProject.filt <- addIterativeLSI(
  ArchRProj = ArchRProject.filt,
  useMatrix = args$matrix, 
  name = "IterativeLSI", 
  firstSelection = "Top",
  depthCol = "nFrags",
  iterations = opts$lsi.iterations, 
  # clusterParams = list(
  #  resolution = opts$lsi.cluster.resolution, 
  #  sampleCells = 10000, 
  #  n.start = 10
  # ), 
  saveIterations = FALSE,
  varFeatures = args$nfeatures, 
  force = TRUE
)

# Correlation between latent dimensions and number of peaks
# cor(getReducedDims(ArchRProject.filt, "IterativeLSI"),ArchRProject.filt$nFrags)[,1] %>% abs %>% sort(decreasing = T)

############################
## LSI + Batch correction ##
############################

if (length(args$batch.variable)>0) {
  print(sprintf("Applying %s batch correction for variable: %s", args$batch.method, args$batch.variable))
  outfile <- sprintf("%s/lsi_%s_nfeatures%d_dims%d_%sbatchcorrection_by_%s.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims, args$batch.method, paste(args$batch.variable,collapse="-"))
  
  # Harmony
  if (args$batch.method=="Harmony") {
    ArchRProject.filt <- addHarmony(
      ArchRProj = ArchRProject.filt,
      reducedDims = "IterativeLSI",
      name = "IterativeLSI_Harmony",
      groupBy = "stage",
      force = TRUE
    )
    lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI_Harmony") %>% round(3) %>% 
      as.data.table(keep.rownames = T) %>% setnames("rn","cell")

  } else {
    stop("Batch correction method not recognised")
  }

} else {
  outfile <- sprintf("%s/lsi_%s_nfeatures%d_ndims%d.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims)
  lsi.dt <- getReducedDims(ArchRProject.filt, "IterativeLSI") %>% round(3) %>% 
    as.data.table(keep.rownames = T) %>% setnames("rn","cell")
}

# Save LSI coordinates
fwrite(lsi.dt, outfile)

##########
## UMAP ##
##########

# i <- args$n_neighbors[1]
# j <- args$min_dist[1]
for (i in args$n_neighbors) {
  for (j in args$min_dist) {
    
    # Define the latent space to run UMAP on
    if (length(opts$batch.correction)>0) {
      if (args$batch.method=="Harmony") {
        dimred <- "IterativeLSI_Harmony"
      } else if  (args$batch.method=="MNN") {
        dimred <- "IterativeLSI_MNN"
      } else {
        stop("Batch correction method not recognised")
      }
    } else {
      dimred <- "IterativeLSI"
    }
    
    # Run UMAP
    ArchRProject.filt <- addUMAP(
      ArchRProj = ArchRProject.filt, 
      reducedDims = dimred,
      name = "UMAP",
      metric = "cosine",
      nNeighbors = i, 
      minDist = j, 
      seed = args$seed,
      saveModel = FALSE,
      force = TRUE
    )
    
    # Fetch UMAP coordinates
    umap.dt <- getEmbedding(ArchRProject.filt,"UMAP") %>%
      round(2) %>%
      as.data.table(keep.rownames = T) %>%
      setnames(c("cell","umap1","umap2"))
    
    # Save UMAP coordinates
    # outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_neigh%d_dist%s.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims, i, j)
    outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d.txt.gz",args$outdir, args$matrix, args$nfeatures, args$ndims)
    fwrite(umap.dt, outfile)

    # Plot
    to.plot <- umap.dt %>%
      merge(sample_metadata,by="cell")
    
    for (k in args$colour_by) {

      # log10 large numeric values
      if (is.numeric(to.plot[[k]])) {
        if (max(to.plot[[k]],na.rm=T) - min(to.plot[[k]],na.rm=T) > 1000) {
          to.plot[[k]] <- log10(to.plot[[k]]+1)
          to.plot %>% setnames(k,paste0(k,"_log10")); k <- paste0(k,"_log10")
        }
      }
      
      p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", fill=k)) +
        geom_point(size=1.5, shape=21, stroke=0.05) +
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
      outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_%s.pdf",args$outdir, args$matrix, args$nfeatures, args$ndims, k)
      pdf(outfile, width=7, height=5)
      print(p)
      dev.off()
    }
    
  }
}

