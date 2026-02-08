source(here::here("settings.R"))
source(here::here("utils.R"))

library(reticulate)

#####################
## Define settings ##
#####################

# I/O
io$metacells.dir <- file.path(io$basedir,"results/metacells/per_condition")
io$metadata <- file.path(io$basedir,"processed/cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/metacells/per_condition")
io$python <- "/share/rargelaguet/mambaforge/envs/metacells/bin/python" 
dir.create(io$outdir, showWarnings = F); dir.create(file.path(io$outdir,"pdf"), showWarnings = F)

# Options
opts$metacell_min_reads <- 25000
opts$metadata.basic.cols <- c("metacell","sample","condition")
opts$metadata.discrete.cols <- c("cell_type")
opts$conditions <- c("Saline Control", "Bleomycin Fibrosis Model", "Bleomycin Fibrosis plus Nintedanib")

# Reticulate
use_python(io$python, required = TRUE)
sc <- import("scanpy")

###########################
## Load metacell results ##
###########################

cell2metacell.dt <- opts$conditions %>% map(function(i) {
  fread(file.path(io$metacells.dir,sprintf("%s/cell2metacell_assignment.txt.gz",gsub(" ","-",i)))) %>% 
    .[,metacell:=paste(i,metacell,sep="#")] %>% .[,condition:=i]
}) %>% rbindlist

stopifnot(unique(cell2metacell.dt$condition)%in%opts$conditions)

# print("Number of metacells per condition:")
# print(table(cell2metacell.dt$condition))
print(sprintf("Total number of metacells: %s", length(unique(cell2metacell.dt$metacell))))

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(io$metadata) %>%
  merge(cell2metacell.dt, by=c("cell","condition"))

cell_metadata.dt %>% setnames("total_counts","number_reads")

###############################
## Load SingleCellExperiment ##
###############################

sce <- load_SingleCellExperiment(file=io$sce, cells=cell_metadata.dt$cell)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

sce <- sce[!duplicated(rownames(sce))]

######################
## Filter metacells ##
######################

# nreads_per_metacell.dt <- data.table(
#   metacell = colnames(sce_metacells),
#   nreads = colSums(counts(sce_metacells))
# )
nreads_per_metacell.dt <- cell_metadata.dt[,.(nreads=sum(number_reads)),by="metacell"]

# Plot
p <- gghistogram(nreads_per_metacell.dt[nreads<=1e6], x="nreads", y="..density..", bins=100, fill="gray70", color="gray50") +
  geom_vline(xintercept=opts$metacell_min_reads, linetype="dashed") +
  labs(x="Number of reads per metacell", y="Density") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank()
  )

pdf(file.path(io$outdir,"pdf/nreads_metacell_threshold.pdf"), width=8, height=4)
print(p)
dev.off()


# Filter metacells
print(sprintf("Removing metacells that have less than %d reads (%d/%d)", opts$metacell_min_reads, nreads_per_metacell.dt[nreads<=opts$metacell_min_reads,.N], nrow(nreads_per_metacell.dt)))
metacells.to.use <- nreads_per_metacell.dt[nreads>=opts$metacell_min_reads,metacell]
cell2metacell.dt <- cell2metacell.dt[metacell%in%metacells.to.use]

# Filter cells
cell_metadata.dt <- cell_metadata.dt[metacell%in%metacells.to.use]
sce <- sce[,cell_metadata.dt$cell]

###################################
## Aggregate counts per metacell ##
###################################

sce_metacells <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = "metacell",
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

# set assay name
assayNames(sce_metacells) <- "counts"

# incorporate cell2metacell assignment
metadata(sce_metacells)$cell2metacell <- cell2metacell.dt

# add rowData
rowData(sce_metacells) <- rowData(sce)

##############################
## Define metacell metadata ##
##############################

metacell_metadata.dt <- data.table(metacell = colnames(sce_metacells), number_reads = colSums(counts(sce_metacells))) %>%
  merge(unique(cell_metadata.dt[,c("metacell","condition")]), by="metacell") %>% 
  setkey(metacell) %>% .[colnames(sce_metacells)]

# Add metacell metadata as colData
stopifnot(metacell_metadata.dt$metacell == colnames(sce_metacells))
colData(sce_metacells) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

############################################
## Infer discrete variables for metacells ##
############################################

metacells_discrete.dt <- data.table(metacell=unique(cell2metacell.dt$metacell))

for (i in opts$metadata.discrete.cols) {
  tmp <- cell2metacell.dt %>% merge(cell_metadata.dt[,c("cell",i),with=F], by="cell")
  tmp.1 <- tmp %>% .[,.N,by=c("metacell",i)] %>% .[,.SD[which.max(N)],by="metacell"] %>% .[,N:=NULL]
  metacells_discrete.dt <- metacells_discrete.dt %>% merge(tmp.1, by="metacell")
}

metacell_metadata.dt <- metacell_metadata.dt %>% merge(metacells_discrete.dt,by="metacell")

stopifnot(metacell_metadata.dt$metacell == colnames(sce_metacells))
colData(sce_metacells) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

###################################################
## Add precomputed dimensionality reduction: PCA ##
###################################################

# io$dimred.cells <- file.path(io$basedir,"results/dimensionality_reduction/precomputed_pca.txt.gz")
# dimred_cells.dt <- fread(io$dimred.cells)# %>% setnames(c("cell","x","y"))
dimred_cells.dt <- reducedDim(sce,"X_pca") %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")

dimred_metacells.dt <- cell2metacell.dt[,c("cell","metacell")] %>%  
  merge(dimred_cells.dt,by="cell") %>% 
  .[, lapply(.SD, mean), by = "metacell", .SDcols = 3:(ncol(dimred_cells.dt)+1)]

reducedDim(sce_metacells,"X_pca") <- dimred_metacells.dt %>% 
  .[metacell%in%colnames(sce_metacells)] %>% setkey(metacell) %>% .[colnames(sce_metacells)] %>%
  matrix.please %>% round(digits=2)

####################################################
## Add precomputed dimensionality reduction: UMAP ##
####################################################

# io$dimred.cells <- file.path(io$basedir,"results/dimensionality_reduction/precomputed_umap.txt.gz")
# dimred_cells.dt <- fread(io$dimred.cells) %>% setnames(c("cell","x","y"))
dimred_cells.dt <- reducedDim(sce,"X_umap") %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")

dimred_metacells.dt <- cell2metacell.dt[,c("cell","metacell")] %>%  
  merge(dimred_cells.dt,by="cell") %>% 
  .[, lapply(.SD, mean), by = "metacell", .SDcols = 3:(ncol(dimred_cells.dt)+1)]

reducedDim(sce_metacells,"X_umap") <- dimred_metacells.dt %>% 
  .[metacell%in%colnames(sce_metacells)] %>% setkey(metacell) %>% .[colnames(sce_metacells)] %>%
  matrix.please %>% round(digits=2)


####################################
## Plot coverage before and after ##
####################################

to.plot <- rbind(
  metacell_metadata.dt[,c("metacell","number_reads")] %>% .[,class:="metacell"],
  cell_metadata.dt[,c("metacell","number_reads")] %>% .[,class:="cell"]
) %>% .[,log_number_reads:=log10(number_reads)] 

p <- gghistogram(to.plot, x="log_number_reads", fill="class", bins=100) +
  labs(x="Number of RNA reads (log10)") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(file.path(io$outdir,"pdf/nreads_before_vs_after.pdf"), width=7, height=4)
print(p)
dev.off()

#############################################
## Convert SingleCellExperiment to anndata ##
#############################################

adata_metacells <- sc$AnnData(
    X   = t(counts(sce_metacells)),
    obs = as.data.frame(colData(sce_metacells)),
    var = data.frame(gene=rownames(sce_metacells), row.names=rownames(sce_metacells)),
    dtype = "float32"
)

# Add cell type colors
# adata_sce$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata_sce$obs$celltype)))])
# adata_sce$uns["celltype_colors"]

# Add dimensionality reduction
adata_metacells$obsm$update(X_pca = reducedDim(sce_metacells,"X_pca"))
adata_metacells$obsm$update(X_umap = reducedDim(sce_metacells,"X_umap"))

##########
## Save ##
##########

adata_metacells$write_h5ad(file.path(io$outdir,"anndata_metacells.h5ad"), compression="gzip")
fwrite(metacell_metadata.dt, file.path(io$outdir,"metacells_metadata.txt.gz"), sep="\t", quote=F, na="NA")
fwrite(cell2metacell.dt, file.path(io$outdir,"cell2metacell_assignment.txt.gz"), sep="\t", quote=F, na="NA")
saveRDS(sce_metacells, file.path(io$outdir,"SingleCellExperiment_metacells.rds"))
