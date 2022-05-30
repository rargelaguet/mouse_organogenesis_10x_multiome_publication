# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$basedir <- file.path(io$basedir,"test")
io$cell_metadata <- io$metadata
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/all_cells/metacells_metadata.txt.gz")
io$pseudobulk_metadata <- file.path(io$basedir,"results/rna/pseudobulk/celltype/stats.txt")
io$cell_sce <- io$rna.sce
io$metacell_sce <- file.path(io$basedir,"results/rna/metacells/all_cells/SingleCellExperiment_metacells.rds")
io$pseudobulk_sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk_with_replicates.rds")
# io$pseudobulk_sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
io$outdir <- file.path(io$basedir,"results/rna/plot_individual_genes/celltypes"); dir.create(io$outdir, showWarnings = F)

# Options

####################
## Load metadata  ##
####################

pseudobulk_metadata.dt <- fread(io$pseudobulk_metadata)

metacell_metadata.dt <- fread(io$metacell_metadata)

cell_metadata.dt <- fread(io$cell_metadata) %>% 
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype)]

###############################
## Load SingleCellExperiment ##
###############################

cells.sce <- load_SingleCellExperiment(io$cell_sce, cells=cell_metadata.dt$cell, normalise = T)
colData(cells.sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

metacells.sce <- load_SingleCellExperiment(io$metacell_sce, cells=metacell_metadata.dt$metacell, normalise = F)
colData(metacells.sce) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

pseudobulk.sce <- readRDS(io$pseudobulk_sce)
pseudobulk.sce$celltype <- colnames(pseudobulk.sce) %>% strsplit("_rep") %>% map_chr(1)

###########
## Parse ##
###########

# Subset celltypes
cell.celltypes <- cell_metadata.dt[,.N,by=c("celltype")] %>% .[N>=1,celltype]
metacell.celltypes <- metacell_metadata.dt[,.N,by=c("celltype")] %>% .[N>=1,celltype]
celltypes.to.plot <- opts$celltypes[opts$celltypes%in%intersect(cell.celltypes,metacell.celltypes)]

# Subset SingleCellExperiment
cells.sce <- cells.sce[,cells.sce$celltype%in%celltypes.to.plot]
metacells.sce <- metacells.sce[,metacells.sce$celltype%in%celltypes.to.plot]
pseudobulk.sce <- pseudobulk.sce[,pseudobulk.sce$celltype%in%celltypes.to.plot]

###########################
## Plot individual genes ##
###########################

give.n <- function(x) { return(c(y = max(x), label = length(x)))}

# genes.to.plot <- c("T","Hoxa9","Tfap2a","Mixl1","Mesp2","Mesp1")
genes.to.plot <- c("Dlx2","Atp10b","Chst3")

# i <- "T"
for (i in genes.to.plot) {
  
  # Plot single cells 
  cells_expr.dt <- data.table(
    cell = colnames(cells.sce),
    expr = logcounts(cells.sce)[i,],
    celltype = cells.sce$celltype
  ) %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  p.cells <- ggplot(cells_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
    # geom_point(shape=21, size=2, data=to.plot[,.(expr=mean(expr)),by="celltype"]) +
    geom_violin(scale="width", alpha=0.75) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75) +
    stat_summary(fun.data = function(x) { return(c(y = max(cells_expr.dt$expr)+0.5, label = length(x)))}, geom = "text", size=2) +
    scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
    labs(x="",y="RNA expression (cells)", title=i) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      # axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  # Plot metacells
  metacells_expr.dt <- data.table(
    cell = colnames(metacells.sce),
    expr = logcounts(metacells.sce)[i,],
    celltype = metacells.sce$celltype
  ) %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  p.metacells <- ggplot(metacells_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
    geom_violin(scale="width", alpha=0.75) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75) +
    stat_summary(fun.data = function(x) { return(c(y = max(metacells_expr.dt$expr)+1, label = length(x)))}, geom = "text", size=2.75) +
    scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
    labs(x="",y="RNA expression (metacells)") +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      # axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  # Plot pseudobulk expr
  pseudobulk_expr.dt <- data.table(
    expr = logcounts(pseudobulk.sce)[i,],
    celltype = pseudobulk.sce$celltype
  ) %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  p.pseudobulk <- ggplot(pseudobulk_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
    geom_bar(stat="identity", color="black", data=pseudobulk_expr.dt[,.(expr=mean(expr)),by="celltype"], alpha=0.75) +
    geom_point(shape=21, size=2) +
    stat_summary(fun.data = function(x) { return(c(y = max(metacells_expr.dt$expr)+1, label = length(x)))}, geom = "text", size=2.75) +
    scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
    labs(x="",y="RNA expression (pseudobulk)") +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      # axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  p <- cowplot::plot_grid(plotlist=list(p.cells,p.metacells, p.pseudobulk), ncol=1, rel_heights = c(0.3,0.3,0.3))
  
  pdf(file.path(io$outdir,sprintf("%s_cell_metacell_pseudobulk_expr.pdf",i)), width=10, height=10)
  print(p)
  dev.off()
}


##########
## TEST ##
##########

# pseudobulk_expr.dt <- data.table(
#   expr = logcounts(pseudobulk.sce)[i,],
#   celltype = colnames(pseudobulk.sce),
#   class = "cell"
# ) %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
# 
# p.pseudobulk <- ggplot(pseudobulk_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
#   geom_bar(stat="identity", color="black") +
#   scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
#   labs(x="",y=sprintf("%s expression (pseudobulk)",i)) +
#   guides(x = guide_axis(angle = 90)) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(colour="black",size=rel(0.8)),
#     axis.text.y = element_text(colour="black",size=rel(0.9)),
#     axis.ticks.x = element_blank(),
#     legend.position = "none"
#   )