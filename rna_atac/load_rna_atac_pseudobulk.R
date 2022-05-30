## START TEST ##
# io$rna.pseudobulk.sce <- file.path(io$basedir,"results_new/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped_mnn.rds")
# io$archR.pseudobulk.deviations.se <- file.path(io$basedir,"results_new/atac/archR/chromvar/pseudobulk/chromVAR_deviations_CISBP_pseudobulk_archr.rds")
# io$archR.peak.metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
# io$archR.peakMatrix.pseudobulk <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn/pseudobulk_PeakMatrix_summarized_experiment.rds")
# io$archR.pseudobulk.GeneScoreMatrix.se <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn/pseudobulk_GeneScoreMatrix_TSS_summarized_experiment.rds")
# 
# opts$motif_annotation <- "CISBP"
## END TEST ##

###################
## Sanity checks ##
###################

if (is.null(opts$motif_annotation)) {
  opts$motif_annotation <- "CISBP"
}

stopifnot(!is.null(opts$celltypes))

stopifnot(file.exists(io$rna.pseudobulk.sce))
stopifnot(file.exists(io$rna.tfs.pseudobulk.sce))
# stopifnot(file.exists(io$archR.pseudobulk.GeneScoreMatrix.se))
stopifnot(file.exists(io$archR.peakMatrix.pseudobulk))
stopifnot(file.exists(io$archR.peak.metadata))
stopifnot(file.exists(io$archR.chromvar.pseudobulk))

##############################################
## Load pseudobulk RNA expression estimates ##
##############################################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]
rna_tf_pseudobulk.sce <- readRDS(io$rna.tfs.pseudobulk.sce)[,opts$celltypes]

##########################################
## Load pseudobulk ATAC GeneScoreMatrix ##
##########################################

# atac_pseudobulk_GeneScoreMatrix.se <- readRDS(io$archR.pseudobulk.GeneScoreMatrix.se)

# Rename genes
# if (any(grepl("f+",rownames(atac_pseudobulk_GeneScoreMatrix.se)))) {
#   rownames(atac_pseudobulk_GeneScoreMatrix.se) <- rowData(atac_pseudobulk_GeneScoreMatrix.se)$name
# }

#####################################
## Load pseudobulk ATAC peakMatrix ##
#####################################

# Load SummarizedExperiment
atac_pseudobulk_peakMatrix.se <- readRDS(io$archR.peakMatrix.pseudobulk)[,opts$celltypes]

# Load peak metadata
# peak_metadata.dt <- fread(io$archR.peak.metadata) %>%
#   .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

# Normalise ATAC data
assayNames(atac_pseudobulk_peakMatrix.se) <- "counts"
assay(atac_pseudobulk_peakMatrix.se,"logcounts") <- log(1e6*(sweep(assay(atac_pseudobulk_peakMatrix.se),2,colSums(assay(atac_pseudobulk_peakMatrix.se),na.rm=T),"/"))+1)

#####################################
## Load pseudobulk chromVAR scores ##
#####################################

# Precomputed
atac_chromvar_pseudobulk.se <- readRDS(io$archR.chromvar.pseudobulk)[,opts$celltypes]

# Set z-scores to the main assay
if ("z"%in%assayNames(atac_chromvar_pseudobulk.se)) {
  assays(atac_chromvar_pseudobulk.se) <- assays(atac_chromvar_pseudobulk.se)["z"]
} else {
  atac_chromvar_pseudobulk.se <- atac_chromvar_pseudobulk.se[rowData(atac_chromvar_pseudobulk.se)$seqnames=="z",]
}


################################
## Load motif2gene annotation ##
################################

motif2gene.dt <- fread(io$archR.motif2gene)

# motif2gene.dt <- motif2gene.dt %>%
#   .[motif%in%rownames(atac_chromvar_pseudobulk.se) & gene%in%toupper(rownames(rna_pseudobulk.sce))] %>%
#   .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

########################
## create data.tables ##
########################

atac_chromvar_pseudobulk.dt <- assay(atac_chromvar_pseudobulk.se) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="motif", value.name="zscore") %>%
  merge(motif2gene.dt[,c("gene","motif")], by="motif", all.x=T)
  
atac_peaks_pseudobulk.dt <- assay(atac_pseudobulk_peakMatrix.se,"logcounts") %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="peak", value.name="acc")

# atac_gene_scores_pseudobulk.dt <- as.matrix(assay(atac_pseudobulk_GeneScoreMatrix.se)) %>% t %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","celltype") %>%
#   melt(id.vars=c("celltype"), variable.name="gene", value.name="acc")

rna_pseudobulk.dt <- logcounts(rna_pseudobulk.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

rna_tf_pseudobulk.dt <- logcounts(rna_tf_pseudobulk.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  data.table::melt(id.vars="gene", variable.name="celltype", value.name="expr")


###################################
## Load chromVAR-ChIP (optional) ##
###################################

if (file.exists(io$archR.chromvar_chip.pseudobulk)) {
  atac_chromvar_chip_pseudobulk.se <- readRDS(io$archR.chromvar_chip.pseudobulk)[,opts$celltypes]
  
  # Set z-scores to the main assay
  if ("z"%in%assayNames(atac_chromvar_chip_pseudobulk.se)) {
    assays(atac_chromvar_chip_pseudobulk.se) <- assays(atac_chromvar_chip_pseudobulk.se)["z"]
  } else {
    atac_chromvar_chip_pseudobulk.se <- atac_chromvar_chip_pseudobulk.se[rowData(atac_chromvar_chip_pseudobulk.se)$seqnames=="z",]
  }
  
  atac_chromvar_chip_pseudobulk.dt <- assay(atac_chromvar_chip_pseudobulk.se) %>% t %>%
    as.data.table(keep.rownames = T) %>%
    setnames("rn","celltype") %>%
    melt(id.vars=c("celltype"), variable.name="gene", value.name="zscore")
}

