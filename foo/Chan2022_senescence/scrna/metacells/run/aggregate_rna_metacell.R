# here::i_am("metacells/run/aggregate_rna_metacell.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

library(reticulate)

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- file.path(io$basedir,"cell_metadata.txt.gz"); io$sce <- file.path(io$basedir,"processed/scTE/TE_genes_SingleCellExperiment.rds"); io$outdir <- file.path(io$basedir,"results/metacells/scTE")
# io$metadata <- file.path(io$basedir,"processed/rna/cellranger/cell_metadata.txt.gz"); io$sce <- file.path(io$basedir,"processed/rna/cellranger/SingleCellExperiment.rds"); io$outdir <- file.path(io$basedir,"results/rna/metacells/cellranger")
io$metadata <- file.path(io$basedir,"results/rna/mapping/cell_metadata_after_mapping.txt.gz"); io$sce <- file.path(io$basedir,"processed/rna/scTE/SingleCellExperiment.rds"); io$outdir <- file.path(io$basedir,"results/rna/metacells/scTE")
io$python <- "/share/rargelaguet/mambaforge/envs/metacells/bin/python"
dir.create(io$outdir, showWarnings = F, recursive=TRUE)
dir.create(file.path(io$outdir,"pdf"), showWarnings = F)

# Options
opts$metacell_min_reads <- 1e4

# Reticulate
use_python(io$python, required = TRUE)
sc <- import("scanpy")

###########################
## Load metacell results ##
###########################

cell2metacell.dt <- opts$samples %>% map(function(i) {
  fread(file.path(io$basedir,sprintf("results/rna/metacells/scTE/%s/cell2metacell_assignment.txt.gz",i))) %>% 
    .[,metacell:=paste(i,metacell,sep="#")]
}) %>% rbindlist

# cell2metacell.dt <- io$cell2metacell %>% map(~ fread(.)) %>% rbindlist

print(sprintf("Number of metacells: %s", length(unique(cell2metacell.dt$metacell))))

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(io$metadata) %>%
  merge(cell2metacell.dt,"cell") %>% setnames("nFeature_RNA","number_reads")

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

metacell_metadata.dt <- data.table(
  metacell = colnames(sce_metacells), 
  number_reads = colSums(counts(sce_metacells))
) %>% merge(unique(cell_metadata.dt[,c("metacell","sample","alias","stage")]), by="metacell")

# Add metacell metadata as colData
stopifnot(metacell_metadata.dt$metacell == colnames(sce_metacells))
colData(sce_metacells) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

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

############################################
## Infer discrete variables for metacells ##
############################################

metacells_discrete.dt <- data.table(metacell=unique(cell2metacell.dt$metacell))

discrete.cols <- c("cell_cycle_phase", "sample_nair2023", "sample_liu2020", "sample_bj_champ")
for (i in discrete.cols) {
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

io$dimred.cells <- file.path(io$basedir,"results/rna/dimensionality_reduction/precomputed_pca.txt.gz")

dimred_cells.dt <- fread(io$dimred.cells)# %>% setnames(c("cell","x","y"))

dimred_metacells.dt <- cell2metacell.dt %>%  
  merge(dimred_cells.dt,by="cell") %>% 
  .[, lapply(.SD, mean), by = "metacell", .SDcols = 3:(ncol(dimred_cells.dt)+1)]

reducedDim(sce_metacells,"PCA") <- dimred_metacells.dt %>% 
  .[metacell%in%colnames(sce_metacells)] %>% setkey(metacell) %>% .[colnames(sce_metacells)] %>%
  matrix.please %>% round(digits=2)

####################################################
## Add precomputed dimensionality reduction: UMAP ##
####################################################

io$dimred.cells <- file.path(io$basedir,"results/rna/dimensionality_reduction/precomputed_umap.txt.gz")

dimred_cells.dt <- fread(io$dimred.cells) %>% setnames(c("cell","x","y"))

dimred_metacells.dt <- dimred_cells.dt %>%  
  merge(cell2metacell.dt,by="cell") %>%
  .[,.(x=mean(x), y=mean(y)), by="metacell"]

reducedDim(sce_metacells,"UMAP") <- dimred_metacells.dt %>% 
  .[metacell%in%colnames(sce_metacells)] %>% setkey(metacell) %>% .[colnames(sce_metacells)] %>%
  matrix.please %>% round(digits=2)

#######################################################
## Add precomputed dimensionality reduction: 3D UMAP ##
#######################################################

# io$dimred.cells <- file.path(io$basedir,"results/dimensionality_reduction/precomputed_umap_3D.txt.gz")

# dimred_cells.dt <- fread(io$dimred.cells) %>% setnames(c("cell","x","y","z"))

# dimred_metacells.dt <- dimred_cells.dt %>%  
#   merge(cell2metacell.dt,by="cell") %>%
#   .[,.(x=mean(x), y=mean(y), z=mean(z)), by="metacell"]

# reducedDim(sce_metacells,"UMAP_3D") <- dimred_metacells.dt %>% 
#   .[metacell%in%colnames(sce_metacells)] %>% setkey(metacell) %>% .[colnames(sce_metacells)] %>%
#   matrix.please %>% round(digits=2)

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
# adata_metacells$uns$update(celltype_colors = opts$celltype.colors[sort(unique(as.character(adata_metacells$obs$celltype)))])
# adata_metacells$uns["celltype_colors"]

##########
## Save ##
##########

adata_metacells$write_h5ad(file.path(io$outdir,"anndata_metacells.h5ad"), compression="gzip")
fwrite(cell2metacell.dt, file.path(io$outdir,"cell2metacell_assignment.txt.gz"), sep="\t", quote=F, na="NA")
fwrite(metacell_metadata.dt, file.path(io$outdir,"metacells_metadata.txt.gz"), sep="\t", quote=F, na="NA")
saveRDS(sce_metacells, file.path(io$outdir,"SingleCellExperiment_metacells.rds"))
