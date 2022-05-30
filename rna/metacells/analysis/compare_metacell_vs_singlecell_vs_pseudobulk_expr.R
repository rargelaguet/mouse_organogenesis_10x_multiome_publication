# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$cell_metadata <- io$metadata
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/metacells_metadata.txt.gz")
io$cell_sce <- io$rna.sce
io$metacell_sce <- file.path(io$basedir,"results/rna/metacells/SingleCellExperiment_metacells.rds")
io$pseudobulk_sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype.mapped/SingleCellExperiment_pseudobulk.rds")
io$pseudobulk_stats <- file.path(io$basedir,"results/rna/pseudobulk/celltype.mapped/stats.txt")
io$outdir <- file.path(io$basedir,"results/rna/metacells/pdf"); dir.create(io$outdir, showWarnings = F)

dir.create(file.path(io$outdir,"individual_genes"), showWarnings = F)

# Options
# opts$samples <- c("E8.75_rep1")

####################
## Load metadata  ##
####################

metacell_metadata.dt <- fread(io$metacell_metadata)

cell_metadata.dt <- fread(io$cell_metadata) %>% 
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype)]

###############################
## Load SingleCellExperiment ##
###############################

cells.sce <- load_SingleCellExperiment(file=io$cell_sce, cells=cell_metadata.dt$cell, remove_non_expressed_genes = T, normalise = T)
colData(cells.sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

metacells.sce <- load_SingleCellExperiment(file=io$metacell_sce, cells=metacell_metadata.dt$metacell, remove_non_expressed_genes = T, normalise = T)
colData(metacells.sce) <- metacell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

pseudobulk.sce <- readRDS(io$pseudobulk_sce)

###########
## Parse ##
###########

# Subset celltypes
cell.celltypes <- cell_metadata.dt[,.N,by=c("celltype")] %>% .[N>=5,celltype]
metacell.celltypes <- metacell_metadata.dt[,.N,by=c("celltype")] %>% .[N>=5,celltype]
celltypes.to.plot <- opts$celltypes[opts$celltypes%in%intersect(cell.celltypes,metacell.celltypes)]

# Subset SingleCellExperiment
cells.sce <- cells.sce[,cells.sce$celltype%in%celltypes.to.plot]
metacells.sce <- metacells.sce[,metacells.sce$celltype%in%celltypes.to.plot]
pseudobulk.sce <- pseudobulk.sce[,colnames(pseudobulk.sce)%in%celltypes.to.plot]

#########################
## Plot coverage stats ##
#########################

to.plot <- rbindlist(
  list(
    data.table(sample=colnames(cells.sce), nFeature_RNA=colSums(counts(cells.sce)), class="cell"),
    data.table(sample=colnames(metacells.sce), nFeature_RNA=colSums(counts(metacells.sce)), class="metacell"),
    data.table(sample=colnames(pseudobulk.sce), nFeature_RNA=colSums(counts(pseudobulk.sce)), class="pseudobulk")
  )
) %>% .[,log2_nFeature_RNA:=log2(nFeature_RNA)] 


p <- ggboxplot(to.plot, x="class", y="log2_nFeature_RNA", fill="gray70") +
  geom_text(aes(label=N, x=class), y=max(to.plot$log2_nFeature_RNA), data=to.plot[,.N,by="class"]) +
  labs(x="", y="Number of RNA reads (log2)") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(file.path(io$outdir,"boxplots_coverage_cell-vs-metacell-vs-pseudobulk.pdf"), width=6, height=4)
print(p)
dev.off()

###########################
## Plot individual genes ##
###########################

give.n <- function(x) { return(c(y = max(x), label = length(x)))}

genes.to.plot <- c("T","Hoxa9","Tfap2a")

# i <- "T"
for (i in genes.to.plot) {
  
  # Plot single cells 
  cells_expr.dt <- data.table(
    cell = colnames(cells.sce),
    expr = logcounts(cells.sce)[i,],
    celltype = cells.sce$celltype,
    class = "cell"
  ) %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  p.cells <- ggplot(cells_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
    # geom_point(shape=21, size=2, data=to.plot[,.(expr=mean(expr)),by="celltype"]) +
    geom_violin(scale="width", alpha=0.75) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75) +
    stat_summary(fun.data = function(x) { return(c(y = max(cells_expr.dt$expr)+0.5, label = length(x)))}, geom = "text", size=2) +
    scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
    labs(x="",y=sprintf("%s expression (cells)",i)) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
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
    celltype = metacells.sce$celltype,
    class = "cell"
  ) %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  p.metacells <- ggplot(metacells_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
    geom_violin(scale="width", alpha=0.75) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75) +
    stat_summary(fun.data = function(x) { return(c(y = max(metacells_expr.dt$expr)+1, label = length(x)))}, geom = "text", size=2.75) +
    scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
    labs(x="",y=sprintf("%s expression (metacells)",i)) +
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
    celltype = colnames(pseudobulk.sce),
    class = "cell"
  ) %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  p.pseudobulk <- ggplot(pseudobulk_expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
    labs(x="",y=sprintf("%s expression (pseudobulk)",i)) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(colour="black",size=rel(0.8)),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  
  p <- cowplot::plot_grid(plotlist=list(p.cells,p.metacells, p.pseudobulk), ncol=1, rel_heights = c(0.3,0.3,0.4))
  
  pdf(file.path(io$outdir,sprintf("individual_genes/%s_cell_metacell_pseudobulk_expr.pdf",i)), width=10, height=10)
  print(p)
  dev.off()
}

