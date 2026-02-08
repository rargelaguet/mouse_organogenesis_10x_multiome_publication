source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(UCell))
suppressPackageStartupMessages(library(msigdbr))

set.seed(42)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"results/rna/metacells/metacells_metadata.txt.gz")
io$sce <- file.path(io$basedir,"results/rna/metacells/SingleCellExperiment_metacells.rds")
io$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction/metacells"); dir.create(io$outdir, showWarnings = F)
io$outdir.signatures <- paste0(io$basedir,"/results/rna/gene_signatures/metacells"); dir.create(io$outdir.signatures, showWarnings = F, recursive = T)
io$outdir.expr <- paste0(io$basedir,"/results/rna/gene_expression/metacells"); dir.create(io$outdir.expr, showWarnings = F, recursive = T)

# Options
opts$samples <- c( 
  "D0",
  "D2",
  "D8",
  "D12",
  "D16_TRA-1-60_negative",
  "D16_TRA-1-60_positive"
)

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[,sample:=NULL] %>% setnames("alias","sample") %>%
  .[sample%in%opts$samples] %>% .[,sample:=factor(sample,levels=opts$samples)] %>%
  .[,stage:=factor(stage,levels=opts$stages)] %>% droplevels %>%
  setnames("metacell","cell")

table(cell_metadata.dt$sample)

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = T)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

# reducedDimNames(sce)
reducedDim(sce,"UMAP")[,1] <- -reducedDim(sce,"UMAP")[,1]
reducedDim(sce,"UMAP")[,2] <- -reducedDim(sce,"UMAP")[,2]
scater::plotReducedDim(sce, dimred="UMAP", colour_by="sample", point_size=1, point_alpha=1)

##########################################
## Plot dimred color by gene expression ##
##########################################

dimred.to.plot <- "UMAP"

genes.to.plot <- c("COL1A1","NANOG","EPCAM","ZEB1","LIN28A","A4GALT","CD75","ST6GAL1","CD177")
genes.to.plot <- c("HERVH-int","HERVK-int","HERVL-int","LTR5_Hs","LTR5A","LTR5B","LTR7A","LTR7B","MLT2A1")
# genes.to.plot <- c("DUXA", "ZSCAN4", "MBD3L2","TRIM43", "RFPL1","LEUTX")
# genes.to.plot <- marker_genes.dt[group=="re4"& score>=0.85,gene] %>% head(n=25)
# genes.to.plot <- marker_genes.dt[group=="c12"& score>=0.85,gene] %>% head(n=25)
genes.to.plot <- unique(marker_genes_liu2020.dt$gene)

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "GDF6"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce)[i,],
    sample = sce$sample,
    x = reducedDim(sce, dimred.to.plot)[,1],
    y = reducedDim(sce, dimred.to.plot)[,2]
  )
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=2, stroke=0.01, alpha=0.90, raster.dpi=150) +
    # facet_wrap(~sample, nrow=3) +
    scale_fill_gradient(low = "gray80", high = "red") +
    # scale_fill_viridis() +
    theme_classic() +
    labs(x="UMAP-1", y="UMAP-2", title=i) +
    # ggplot_theme_NoAxes() +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust=0.5),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  # pdf(file.path(io$outdir.expr,sprintf("Xing2020_%s_%s_expr_cells.pdf",dimred.to.plot,i)), width=6, height=9)
  pdf(file.path(io$outdir.expr,sprintf("%s_expr_%s_Xing2020_metacells.pdf",i,dimred.to.plot)), width=6, height=6)
  print(p)
  dev.off()
}
