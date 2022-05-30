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
io$rna.sce  <- paste0(io$basedir,"/processed/rna/SingleCellExperiment_velocyto.rds")
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
sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples & !is.na(celltype.mapped)]

# Load velocyto SingleCellExperiment
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$cell)

# Filter genes
# gene_metadata <- fread(io$gene_metadata)
# genes <- unique(gene_metadata$symbol)
# sce <- sce[rownames(sce)%in%genes,]
sce <- sce[!duplicated(rownames(sce)),]

###################################
## Aggregate counts per celltype ##
###################################

sce_pseudobulk_unspliced <- aggregateData(
  sce,
  assay = "unspliced",
  by = c("celltype.mapped"),
  fun = c("sum"),
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
); assayNames(sce_pseudobulk_unspliced) <- "counts"

sce_pseudobulk_spliced <- aggregateData(
  sce,
  assay = "spliced",
  by = c("celltype.mapped"),
  fun = c("sum"),
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
); assayNames(sce_pseudobulk_spliced) <- "counts"


###############
## Normalise ##
###############

# create DESeq object
dds.unspliced <- DESeqDataSet(sce_pseudobulk_unspliced, design=~1)
dds.spliced <- DESeqDataSet(sce_pseudobulk_spliced, design=~1)

# This function calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) 
# and then transforms the count data (normalized by division by the size factors or normalization factors), 
# yielding a matrix of values which are now approximately homoskedastic 
dds.unspliced <- varianceStabilizingTransformation(dds.unspliced)
dds.spliced <- varianceStabilizingTransformation(dds.spliced)

sce_pseudobulk <- SingleCellExperiment(
  assays = list(
    "spliced" = assay(sce_pseudobulk_spliced), 
    "unspliced" = assay(sce_pseudobulk_unspliced),
    "unspliced_log" = assay(dds.unspliced),
    "spliced_log" = assay(dds.spliced)
    )
)

colnames(sce_pseudobulk) %>% head
rownames(sce_pseudobulk) %>% head

##########
## Save ##
##########

saveRDS(sce_pseudobulk, paste0(io$outdir,"/SingleCellExperiment_velocyto.rds"))
