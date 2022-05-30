here::i_am("rna/metacells/run/aggregate_rna_metacell.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

library(reticulate)

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--sce',    type="character",  help='SingleCellExperiment')
p$add_argument('--python',   type="character",    help='Python path for reticulate')
p$add_argument('--cell2metacell',    type="character",  nargs="+", help='Metacell results')
p$add_argument('--metacell_min_reads',     type="integer",    default=1e5,     help='Minimum number of reads per metacell')
p$add_argument('--normalisation_method',    type="character",    help='Metadata column to group cells by')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
# # args$cell2metacell <- file.path(io$basedir,sprintf("results/rna/metacells/%s/cell2metacell_assignment.txt.gz",opts$samples))
# # args$cell2metacell <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/cell2metacell_assignment.txt.gz")
# args$cell2metacell <- file.path(io$basedir,"results/rna/metacells/all_cells/E8.75_rep1/cell2metacell_assignment.txt.gz")
# args$python <- "/bi/group/reik/ricard/software/miniconda3/envs/main/bin/python" # "/Users/argelagr/opt/anaconda3/envs/main/bin/python"
# args$metacell_min_reads <- 2e4
# args$normalisation_method <- "cpm"
# args$outdir <- file.path(io$basedir,"results/rna/metacells/test")
## END TEST ##

#####################
## Define settings ##
#####################

stopifnot(file.exists(args$cell2metacell))
dir.create(args$outdir, showWarnings = F)

# Reticulate
use_python(args$python, required = TRUE)
sc <- import("scanpy")

###########################
## Load metacell results ##
###########################

cell2metacell.dt <- args$cell2metacell %>% map(~ fread(.)) %>% rbindlist
# stopifnot(mapping_mnn.dt$cell%in%sample_metadata$cell)

print(sprintf("Number of metacells: %s", length(unique(cell2metacell.dt$metacell))))

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>%
  # .[cell%in%cell2metacell.dt$cell] %>%
  # .[,c("cell","sample","stage","genotype","celltype.mapped","closest.cell")] %>%
  # setnames("celltype.mapped","celltype") %>%
  merge(cell2metacell.dt,"cell")

###############################
## Load SingleCellExperiment ##
###############################

sce <- load_SingleCellExperiment(file=args$sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

table(sce$metacell)

######################
## Filter metacells ##
######################

# nreads_per_metacell.dt <- data.table(
#   metacell = colnames(sce_pseudobulk),
#   nreads = colSums(counts(sce_pseudobulk))
# )
nreads_per_metacell.dt <- sample_metadata[,.(nreads=sum(nFeature_RNA)),by="metacell"]

# Plot
p <- gghistogram(nreads_per_metacell.dt[nreads<=1e6], x="nreads", y="..density..", bins=100, fill="gray70", color="gray50") +
  geom_vline(xintercept=args$metacell_min_reads, linetype="dashed") +
  labs(x="Number of reads per metacell", y="Density") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank()
  )

pdf(file.path(args$outdir,"nreads_metacell_threshold.pdf"), width=8, height=4)
print(p)
dev.off()


# Filter
print(sprintf("Removing metacells that have less than %d reads (%d/%d)", args$metacell_min_reads, nreads_per_metacell.dt[nreads<=args$metacell_min_reads,.N], nrow(nreads_per_metacell.dt)))

metacells.to.use <- intersect(sample_metadata$cell, nreads_per_metacell.dt[nreads>=args$metacell_min_reads,metacell])
sample_metadata <- sample_metadata[metacell%in%metacells.to.use]
cell2metacell.dt <- cell2metacell.dt[metacell%in%metacells.to.use]
sce <- sce[,colnames(sce)%in%sample_metadata$cell]

###################################
## Aggregate counts per metacell ##
###################################

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = "metacell",
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

# set assay name
assayNames(sce_pseudobulk) <- "counts"

###################
## Normalisation ##
###################

if (args$normalisation_method=="deseq2") {
  
  suppressPackageStartupMessages(library(DESeq2))
  dds <- DESeqDataSet(sce_pseudobulk, design=~1)
  dds <- varianceStabilizingTransformation(dds)
  logcounts(sce_pseudobulk) <- assay(dds)
  
} else if (args$normalisation_method=="cpm") {
  
  logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)
  # atac.mtx <- log2(t(t(assay(atac.se)) / colSums(assay(atac.se))) * 1e6)
  
} else {
  stop("Normalisation method not recognised")
}

##############################
## Define metacell metadata ##
##############################

metacell_metadata.dt <- sample_metadata %>%
  .[,c("cell","sample","stage","genotype","celltype","closest.cell")] %>%
  .[cell%in%colnames(sce_pseudobulk)] %>% setkey(cell) %>% .[colnames(sce_pseudobulk)] %>%
  setnames("cell","metacell") %>%
  merge(cell2metacell.dt[,.(ncells=.N),by="metacell"],by="metacell")

# Add QC stats
tmp <- data.table(
  metacell = colnames(sce_pseudobulk), 
  nFeature_RNA = colSums(counts(sce_pseudobulk))
)
metacell_metadata.dt <- metacell_metadata.dt %>% merge(tmp,by="metacell")

# Add metacell metadata as colData
stopifnot(metacell_metadata.dt$metacell == colnames(sce_pseudobulk))
colData(sce_pseudobulk) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

####################################
## Plot coverage before and after ##
####################################

# to.plot <- rbind(
#   metacell_metadata.dt[,c("metacell","nFeature_RNA")] %>% .[,class:="metacell"],
#   sample_metadata[,c("metacell","nFeature_RNA")] %>% .[,class:="cell"]
# ) %>% .[,log_nFeature_RNA:=log10(nFeature_RNA)] 

# p <- gghistogram(to.plot, x="log_nFeature_RNA", fill="class", bins=100) +
#   labs(x="Number of RNA reads (log10)") +
#   theme(
#     legend.title = element_blank(),
#     axis.text = element_text(size=rel(0.75), color="black")
#   )

# pdf(file.path(args$outdir,"nreads_before_vs_after.pdf"), width=7, height=4)
# print(p)
# dev.off()

#############################################
## Convert SingleCellExperiment to anndata ##
#############################################

adata_sce <- sc$AnnData(
    X   = t(counts(sce_pseudobulk)),
    obs = as.data.frame(colData(sce_pseudobulk)),
    var = data.frame(gene=rownames(sce_pseudobulk), row.names=rownames(sce_pseudobulk))
)


# Add cell type colors
# adata_sce$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata_sce$obs$celltype)))])
# adata_sce$uns["celltype_colors"]

##########
## Save ##
##########

adata_sce$write_h5ad(file.path(args$outdir,"anndata_metacells.h5ad"))
fwrite(metacell_metadata.dt, file.path(args$outdir,"metacells_metadata.txt.gz"), sep="\t", quote=F, na="NA")
saveRDS(sce_pseudobulk, file.path(args$outdir,"SingleCellExperiment_metacells.rds"))
