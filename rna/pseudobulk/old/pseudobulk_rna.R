library(muscat)
library(DESeq2)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna/pseudobulk")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

###############
## Load data ##
###############

# Load cell metadata
# io$metadata <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/doublets/sample_metadata_after_doublets.txt.gz"
sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples & !is.na(celltype.mapped)]

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

###################################
## Aggregate counts per celltype ##
###################################

# assays(sce)$cpm <- edgeR::cpm(assay(sce), normalized.lib.sizes = FALSE, log = FALSE)

sce_pseudobulk <- aggregateData(
  sce,
  assay = "counts",
  by = c("celltype.mapped"),
  fun = c("sum"),
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

###############
## Normalise ##
###############

# create DESeq object
dds <- DESeqDataSet(sce_pseudobulk, design=~1)

# This function calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) 
# and then transforms the count data (normalized by division by the size factors or normalization factors), 
# yielding a matrix of values which are now approximately homoskedastic 
dds <- varianceStabilizingTransformation(dds)

logcounts(sce_pseudobulk) <- assay(dds)

###################
## Sanity checks ##
###################

# cor(
#   colMeans(logcounts(sce_pseudobulk)),
#   metadata(sce_pseudobulk)$n_cells
# )

##########
## Save ##
##########

saveRDS(sce_pseudobulk, paste0(io$outdir,"/SingleCellExperiment.rds"))
