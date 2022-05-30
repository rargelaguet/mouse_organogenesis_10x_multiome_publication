# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$cell_metadata <- io$metadata
io$cell_sce <- io$rna.sce
io$cell_umap <- file.path(io$basedir,"results/rna/dimensionality_reduction/sce/batch_correction_by_sample_remove_ExE_cells_False/umap_features2500_pcs50_neigh25_dist0.5.txt.gz")
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/metacells_metadata.txt.gz")
io$metacell_sce <- file.path(io$basedir,"results/rna/metacells/SingleCellExperiment_metacells.rds")
io$metacell_umap <- file.path(io$basedir,"results/rna/dimensionality_reduction/metacells/umap_features2500_pcs50_neigh25_dist0.5.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/metacells/individual_genes"); dir.create(io$outdir, showWarnings = F)

###################
## Load metadata ##
###################

metacell_metadata.dt <- fread(io$metacell_metadata) %>% .[,cell:=NULL]

cell_metadata.dt <- fread(io$cell_metadata) %>% 
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype)]

################
## Load UMAPs ##
################

cell_umap.dt <- fread(io$cell_umap, select=c(3,1,2)) %>% setnames(c("cell","V1","V2"))
metacell_umap.dt <- fread(io$metacell_umap, select=c(3,1,2)) %>% setnames(c("metacell","V1","V2"))

# Filter common cells
cells.to.use <- intersect(cell_metadata.dt$cell,cell_umap.dt$cell)
cell_umap.dt <- cell_umap.dt[cell%in%cells.to.use]
cell_metadata.dt <- cell_metadata.dt[cell%in%cells.to.use]

# Filter common metacells
metacells.to.use <- intersect(metacell_metadata.dt$cell,metacell_umap.dt$cell)
metacell_umap.dt <- metacell_umap.dt[cell%in%metacells.to.use]
metacell_metadata.dt <- metacell_metadata.dt[metacell%in%metacells.to.use]

###############################
## Load SingleCellExperiment ##
###############################

cells.sce <- load_SingleCellExperiment(file=io$cell_sce, cells=cell_metadata.dt$cell, remove_non_expressed_genes = T, normalise = T)
colData(cells.sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

metacells.sce <- readRDS(io$metacell_sce)[,metacell_metadata.dt$metacell]
colData(metacells.sce) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

##########
## Plot ##
##########

genes.to.plot <- c("T","Hoxa9","Tfap2a")

# i <- "T"
for (i in genes.to.plot) {
  
  # Plot single cells 
  cells_expr.dt <- data.table(
    cell = colnames(cells.sce),
    celltype = cells.sce$celltype,
    expr = logcounts(cells.sce)[i,]
  ) %>% merge(cell_umap.dt,by="cell")
  
  cells_expr.dt <- cells_expr.dt[sample(1:.N, size=.N/3)]
  
  p.cells.celltype <- ggplot(cells_expr.dt, aes(x=V1, y=V2, color=celltype)) +
    scale_color_manual(values=opts$celltype.colors) +
    ggrastr::geom_point_rast(size=0.25) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position = "none"
    )
  
  p.cells.expr <- ggplot(cells_expr.dt, aes(x=V1, y=V2, color=expr)) +
    scale_color_gradient(low = "gray80", high = "red") +
    geom_point(size=0.25) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position = "none"
    )
  
  # Plot metacells 
  metacells_expr.dt <- data.table(
    metacell = colnames(metacells.sce),
    celltype = metacells.sce$celltype,
    expr = logcounts(metacells.sce)[i,]
  ) %>% merge(metacell_umap.dt,by="metacell")
  
  p.metacells.celltype <- ggplot(metacells_expr.dt, aes(x=V1, y=V2, color=celltype)) +
    scale_color_manual(values=opts$celltype.colors) +
    ggrastr::geom_point_rast(size=0.5) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position = "none"
    )
  
  p.metacells.expr <- ggplot(metacells_expr.dt, aes(x=V1, y=V2, color=expr)) +
    scale_color_gradient(low = "gray80", high = "red") +
    geom_point(size=0.5) +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position = "none"
    )
  
  
  p <- cowplot::plot_grid(plotlist=list(p.cells.celltype,p.cells.expr,p.metacells.celltype,p.metacells.expr), ncol=2)
  
  pdf(file.path(io$outdir,sprintf("%s_umap_cell_metacell_expr.pdf",i)), width=10, height=8)
  print(p)
  dev.off()
}

