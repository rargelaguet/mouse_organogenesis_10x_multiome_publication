here::i_am("rna/metacells/run/aggregate_rna_metacell.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

library(reticulate)
library(zellkonverter)

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--sce',    type="character",  help='SingleCellExperiment')
p$add_argument('--python',   type="character",    help='Python path for reticulate')
p$add_argument('--cell2metacell',    type="character",  nargs="+", help='Metacell results')
p$add_argument('--metacell_min_reads',     type="integer",    default=1e5,     help='Minimum number of reads per metacell')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$sce <- file.path(io$basedir,"processed/rna/velocyto/SingleCellExperiment_velocyto.rds")
args$cell2metacell <- file.path(io$basedir,sprintf("results/rna/metacells/all_cells/%s/cell2metacell_assignment.txt.gz",opts$samples))
# args$cell2metacell <- file.path(io$basedir,"results/rna/metacells/all_cells/E8.75_rep1/cell2metacell_assignment.txt.gz")
args$python <- "/bi/group/reik/ricard/software/miniconda3/envs/main/bin/python" # "/Users/argelagr/opt/anaconda3/envs/main/bin/python"
args$metacell_min_reads <- 2.5e4
args$outdir <- file.path(io$basedir,"results/rna/velocyto/metacells")
## END TEST ##

#####################
## Define settings ##
#####################

stopifnot(file.exists(args$cell2metacell))
dir.create(args$outdir, showWarnings = F, recursive = T)

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

print("Loading SingleCellExperiment object...")

sce <- load_SingleCellExperiment(file=args$sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

table(sce$metacell)

######################
## Filter metacells ##
######################

print("Filtering metacells...")

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

print("Aggregating counts from cells to metacells...")

counts_pseudobulk.sce <- pseudobulk_sce_fn(sce, assay = "counts", by = "metacell", fun = "sum")
spliced_pseudobulk.sce <- pseudobulk_sce_fn(sce, assay = "spliced", by = "metacell", fun = "sum")
unspliced_pseudobulk.sce <- pseudobulk_sce_fn(sce, assay = "unspliced", by = "metacell", fun = "sum")

# set assay name
assayNames(counts_pseudobulk.sce) <- "counts"
assay(counts_pseudobulk.sce, "spliced") <- assay(spliced_pseudobulk.sce)
assay(counts_pseudobulk.sce, "unspliced") <- assay(unspliced_pseudobulk.sce)

##############################
## Define metacell metadata ##
##############################

metacell_metadata.dt <- sample_metadata %>%
  .[,c("cell","sample","stage","genotype","celltype","closest.cell")] %>%
  .[cell%in%colnames(counts_pseudobulk.sce)] %>% setkey(cell) %>% .[colnames(counts_pseudobulk.sce)] %>%
  setnames("cell","metacell") %>%
  merge(cell2metacell.dt[,.(ncells=.N),by="metacell"],by="metacell")

# Add QC stats
tmp <- data.table(
  metacell = colnames(counts_pseudobulk.sce), 
  nFeature_RNA = colSums(counts(counts_pseudobulk.sce))
)
metacell_metadata.dt <- metacell_metadata.dt %>% merge(tmp,by="metacell")

# Add metacell metadata as colData
stopifnot(metacell_metadata.dt$metacell == colnames(counts_pseudobulk.sce))
colData(counts_pseudobulk.sce) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

#############################################
## Convert SingleCellExperiment to anndata ##
#############################################

print("Converting SingleCellExperiment to anndata...")

# data <- sc$AnnData(
#     X   = t(counts(counts_pseudobulk.sce)),
#     # layers = list("spliced"=,"")
#     obs = as.data.frame(colData(counts_pseudobulk.sce)),
#     var = data.frame(gene=rownames(counts_pseudobulk.sce), row.names=rownames(counts_pseudobulk.sce))
# )

adata <- SCE2AnnData(counts_pseudobulk.sce)

# Add cell type colors
adata$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata$obs$celltype)))])

###################
## Normalisation ##
###################

logcounts(counts_pseudobulk.sce) <- log2(1e6*(sweep(counts(counts_pseudobulk.sce),2,colSums(counts(counts_pseudobulk.sce)),"/"))+1)

##########
## Save ##
##########

adata$write_h5ad(file.path(args$outdir,"anndata_metacells.h5ad"), compression="gzip")
fwrite(metacell_metadata.dt, file.path(args$outdir,"metacells_metadata.txt.gz"), sep="\t", quote=F, na="NA")
saveRDS(counts_pseudobulk.sce, file.path(args$outdir,"SingleCellExperiment_metacells.rds"))
