
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/gene_expression"); dir.create(io$outdir, showWarnings = F)

# Options
opts$conditions <- c("Saline_Control", "Bleomycin_day11", "Bleomycin_nintedanib_day11")

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] %>% .[,condition:=factor(condition,levels=opts$conditions)] %>%
  .[cell_type%in%opts$celltypes] %>% .[,cell_type:=factor(cell_type,levels=opts$celltypes)]

table(cell_metadata.dt$condition)

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = T)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

####################
## Load gene sets ##
####################

# gene_metadata.dt <- fread(io$gene_metadata)

##########
## UMAP ##
##########

genes.to.plot <- c("Pi16","Cd34","Has2","Has2","Klf4","Cd70","St6gal1","Cd24a","Podxl","Anpep")
genes.to.plot <- c("St6gal1","Cd24a","Podxl","Anpep","Cd177")
genes.to.plot <- c("Pou5f1","Klf4","Sox2","Myc")
genes.to.plot <- grep("Ly6",rownames(sce),value=T)

genes.to.plot <- fread(io$TFs)[[1]] %>% str_to_title()

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

genes.to.plot <- c("Ly6a","Pi16")

# i <- "S100a8"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = assay(sce,"logcounts")[i,],
    x = reducedDim(sce, "X_umap")[,1],
    y = reducedDim(sce, "X_umap")[,2]
  ) %>% merge(cell_metadata.dt,by="cell")
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=0.75, stroke=0.01, alpha=0.90, raster.dpi=150) +
    # scale_fill_gradient(low = "gray90", high = "red") +
    scale_fill_viridis() +
    theme_classic() +
    labs(x="UMAP-1", y="UMAP-2", title=i) +
    # ggplot_theme_NoAxes() +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust=0.5, size=rel(2)),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  pdf(file.path(io$outdir,sprintf("%s_expr_umap_Peyser2019_fibrosis_cells.pdf",i)), width=5.5, height=5.5)
  print(p)
  dev.off()
}

##############################
## UMAP, split by condition ##
##############################

# i <- "Fcgrt"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = assay(sce,"logcounts")[i,],
    x = reducedDim(sce, "X_umap")[,1],
    y = reducedDim(sce, "X_umap")[,2]
  ) %>% merge(cell_metadata.dt,by="cell")
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=1, stroke=0.01, alpha=0.90, raster.dpi=150) +
    facet_wrap(~condition, nrow=1) +
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
      strip.text = element_text(size=rel(1.10)),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  # pdf(file.path(io$outdir,sprintf("%s_expr_boxplots_Peyser2019_fibrosis_cells.pdf",i)), width=9.5, height=6.5)
  print(p)
  # dev.off()
}


##############################
## UMAP, split by sample ##
##############################

# i <- "Mt1"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = assay(sce,"logcounts")[i,],
    x = reducedDim(sce, "X_umap")[,1],
    y = reducedDim(sce, "X_umap")[,2]
  ) %>% merge(cell_metadata.dt,by="cell")
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=0.75, stroke=0.01, alpha=0.90, raster.dpi=150) +
    facet_wrap(~sample, nrow=3) +
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
      strip.text = element_text(size=rel(0.60)),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  # pdf(file.path(io$outdir,sprintf("%s_expr_boxplots_Peyser2019_fibrosis_cells.pdf",i)), width=9.5, height=6.5)
  print(p)
  # dev.off()
  
}


###############
## Box plots ##
###############

genes.to.plot <- gene_metadata.dt[chr=="chrMT",symbol]

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "Klf4"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = assay(sce,"logcounts")[i,]
  ) %>% merge(cell_metadata.dt,by="cell") 
  
  p <- ggplot(to.plot, aes(x=cell_type, y=expr, fill=cell_type)) +
    # geom_violin(fill="gray80") + geom_boxplot(outlier.shape=NA, fill="gray80") +
    geom_violin(show.legend=F) + geom_boxplot(outlier.shape=NA) +
    stat_summary(fun.data = give.n, fun.args = list(ymax = max(to.plot$expr)+1), size=2.5, geom = "text", position = position_dodge(width=0.75)) +
    facet_wrap(~condition, scales="free_x", nrow=3) +
    # facet_wrap(~donor_id, scales="free_x", space="free_x") +
    guides(x = guide_axis(angle = 90)) +
    labs(x="", y="gene expression (log counts)", title=i) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size=rel(1.10)),
      axis.text.y = element_text(color="black", size=rel(1)),
      # axis.text.x = element_text(color="black", size=rel(1.10)),
      axis.text.x = element_blank(),
      axis.title = element_text(color="black", size=rel(1.10)),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  pdf(file.path(io$outdir,sprintf("%s_expr_boxplots_Peyser2019_fibrosis_cells.pdf",i)), width=9.5, height=6.5)
  print(p)
  dev.off()
}

