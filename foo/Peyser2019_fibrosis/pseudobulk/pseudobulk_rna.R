# here::i_am("rna/pseudobulk/pseudobulk_rna.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$samples <- c("25F", "28F", "52F", "54F", "63F", "70M", "FSFA", "FSFB")
opts$group_by <- "sample"
opts$normalisation_method <- "cpm"
opts$log <- TRUE
opts$min_cells <- 25

# I/O
io$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
io$metadata <- file.path(io$basedir,"processed/rna/cell_metadata.txt.gz")
io$outdir <- file.path(io$basedir,sprintf("results/rna/pseudobulk/%s",opts$group_by)); dir.create(io$outdir, showWarnings = F, recursive = T)

###############
## Load data ##
###############

# Load cell metadata
cell_metadata.dt <- fread(io$metadata) %>% .[sample%in%opts$samples]

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

########################
## Filter by coverage ##
########################

pseudobulk_stats.dt <- data.table(table(cell_metadata.dt[[opts$group_by]])) %>% 
  setnames(c("sample_name","ncells")) 

pseudobulk_stats_no.dt <- pseudobulk_stats.dt %>% .[ncells<opts$min_cells]
pseudobulk_stats_yes.dt <- pseudobulk_stats.dt %>% .[ncells>=opts$min_cells]
cell_metadata.dt <- cell_metadata.dt[cell_metadata.dt[[opts$group_by]]%in%pseudobulk_stats_yes.dt$sample_name]
sce <- sce[,cell_metadata.dt$cell]

# print
print("Groups with sufficient coverage:")
print(pseudobulk_stats_yes.dt)

print("Groups with insufficient coverage:")
print(pseudobulk_stats_no.dt)

################
## Pseudobulk ##
################

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = opts$group_by,
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

###################
## Normalisation ##
###################

if (opts$log) {
  logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)
} else {
  logcounts(sce_pseudobulk) <- 1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))
}

#####################
## Sample metadata ##
#####################

pseudobulk_metadata.dt <- cell_metadata.dt[,c("sample")] %>% unique %>%
  setnames("sample","sample_name") %>%
  merge(pseudobulk_stats.dt[,c("sample_name","ncells")], by="sample_name")

stopifnot(sort(colnames(sce_pseudobulk))==sort(pseudobulk_metadata.dt$sample_name))

colData(sce_pseudobulk) <- pseudobulk_metadata.dt %>% tibble::column_to_rownames("sample_name") %>% DataFrame 

##########
## Save ##
##########

fwrite(pseudobulk_metadata.dt, file.path(io$outdir,"sample_metadata_pseudobulk.txt"), sep="\t", quote = F)
saveRDS(sce_pseudobulk, file.path(io$outdir,"SingleCellExperiment_pseudobulk.rds"))
