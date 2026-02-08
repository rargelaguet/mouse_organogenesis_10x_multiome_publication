source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$group_by <- "alias"
opts$min_cells <- 50

# I/O
io$metadata <- file.path(io$basedir,"cell_metadata.txt.gz")
io$sample_metadata <- file.path(io$basedir,"sample_metadata.txt")
io$sce <- file.path(io$basedir,"processed/cellranger/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/pseudobulk"); dir.create(io$outdir, showWarnings = F, recursive = T)

###############
## Load data ##
###############

# Load cell metadata
cell_metadata.dt <- fread(io$metadata) %>% 
  .[!is.na(eval(as.name(opts$group_by)))]

# Load SingleCellExperiment
sce <- readRDS(io$sce)[,cell_metadata.dt$cell]
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

########################
## Filter by coverage ##
########################

pseudobulk_stats.dt <- data.table(table(cell_metadata.dt[[opts$group_by]])) %>% 
  setnames(c("group","ncells")) 

pseudobulk_stats_no.dt <- pseudobulk_stats.dt %>% .[ncells<opts$min_cells]
pseudobulk_stats_yes.dt <- pseudobulk_stats.dt %>% .[ncells>=opts$min_cells]
cell_metadata.dt <- cell_metadata.dt[cell_metadata.dt[[opts$group_by]]%in%pseudobulk_stats_yes.dt$group]
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
  scale = FALSE
)

assayNames(sce_pseudobulk) <- "counts"

rowData(sce_pseudobulk) <- rowData(sce)

############################
## Update sample metadata ##
############################

sample_metadata.dt <- fread(io$sample_metadata) 

stopifnot(sort(colnames(sce_pseudobulk))==sort(sample_metadata.dt$alias))

colData(sce_pseudobulk) <- sample_metadata.dt %>% setkey(alias) %>%
  .[colnames(sce_pseudobulk)] %>%
  tibble::column_to_rownames("alias") %>% DataFrame

###################
## Normalisation ##
###################

assay(sce_pseudobulk,"cpm") <- 1e6*(sweep(assay(sce_pseudobulk,"counts"),2,colSums(assay(sce_pseudobulk,"counts")),"/"))
assay(sce_pseudobulk,"log_cpm") <- log2(assay(sce_pseudobulk,"cpm")+1)

##########
## Save ##
##########

fwrite(pseudobulk_stats.dt, file.path(io$outdir,"stats.txt"), sep="\t", quote = F)
saveRDS(sce_pseudobulk, file.path(io$outdir,"SingleCellExperiment_pseudobulk.rds"))

# Subset controls
# saveRDS(sce_pseudobulk[,sce_pseudobulk$alias=="Control"], file.path(io$outdir,"SingleCellExperiment_pseudobulk_controls.rds"))

