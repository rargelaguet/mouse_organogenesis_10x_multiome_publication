# here::i_am("dimensionality_reduction/dimensionality_reduction_sce.R")

library(scran)
library(scater)

source(here::here("settings.R"))
source(here::here("utils.R"))

set.seed(42)

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"processed/cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction"); dir.create(io$outdir, showWarnings = F)
io$outdir.expr <- paste0(io$basedir,"/results/gene_expression/fibroblasts"); dir.create(io$outdir.expr, showWarnings = F)

# Options
# opts$conditions <- c("Saline Control", "Bleomycin Fibrosis Model", "Bleomycin Fibrosis plus Nintedanib")
opts$conditions <- c("Saline Control", "Bleomycin Fibrosis Model")

opts$celltypes <- c("fibroblasts")

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] %>% .[,condition:=factor(condition,levels=opts$conditions)] %>%
  .[cell_type%in%opts$celltypes]

cell_metadata.dt[,condition_rename:=as.character(NA)]
cell_metadata.dt[condition=="Saline Control",condition_rename:="Saline_Control"]
cell_metadata.dt[condition=="Bleomycin Fibrosis Model",condition_rename:="Bleomycin_day11"]
cell_metadata.dt[condition=="Bleomycin Fibrosis plus Nintedanib",condition_rename:="Bleomycin_nintedanib_day11"]
cell_metadata.dt[,condition_rename:=factor(condition_rename,levels=c("Saline_Control","Bleomycin_day11","Bleomycin_nintedanib_day11"))]
cell_metadata.dt <- cell_metadata.dt %>% droplevels

table(cell_metadata.dt$condition_rename)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

###############################################
## Load precomputed dimensionality reduction ##
###############################################

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=2500) %>% rownames
sce_filt <- sce[hvgs,]

#########
## PCA ##
#########

sce_filt <- runPCA(sce_filt, ncomponents=30, ntop=2500)
reducedDims(sce) <- reducedDims(sce_filt)

# plotPCA(sce, colour_by="sample", point_size=1, point_alpha=1, ncomponents = c(1,2))
# plotPCA(sce, colour_by="TOP2A", point_size=1, point_alpha=1)

##########
## UMAP ##
##########

sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 25, min_dist = 0.3)
reducedDims(sce) <- reducedDims(sce_filt)
# 
plotReducedDim(sce, dimred="UMAP", colour_by="sample", point_size=1, point_alpha=1)
plotReducedDim(sce, dimred="UMAP", colour_by="condition", point_size=1, point_alpha=1)
plotReducedDim(sce, dimred="UMAP", colour_by="celltype_liu2020", point_size=1, point_alpha=1)
plotReducedDim(sce, dimred="UMAP", colour_by="mit_percent_RNA", point_size=1, point_alpha=1)
plotReducedDim(sce, dimred="UMAP", colour_by="nFeature_RNA", point_size=1, point_alpha=1)
# plotReducedDim(sce, dimred="UMAP", colour_by="HERVK-int", point_size=1, point_alpha=1)
plotReducedDim(sce, dimred="UMAP", colour_by="HERVK-int", point_size=1, point_alpha=1)

######################
## Plot all samples ##
######################

dimred.to.plot <- "PCA"

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(fill=condition), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  ggrepel::geom_text_repel(aes(label=condition), data=to.plot[,.(V1=median(V1), V2=median(V2)), by=c("condition")]) +
  guides(fill = guide_legend(override.aes = list(size=3))) +
  labs(x="UMAP-1", y="UMAP-2") +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  # ggplot_theme_NoAxes() +
  theme(
    strip.text = element_text(size=rel(1.15)),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position  ="none",
    legend.direction = "vertical",
    legend.title=element_blank()
  )

pdf(file.path(io$outdir,sprintf("%s_cell_cycle_Peyser2019_cells.pdf",dimred.to.plot)), width=6, height=5)
print(p)
dev.off()

############################
## Plot, split conditions ##
############################

dimred.to.plot <- "UMAP"

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(fill=cell_type), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  ggrepel::geom_text_repel(aes(label=cell_type), data=to.plot[,.(V1=median(V1), V2=median(V2)), by=c("cell_type")]) +
  facet_wrap(~condition, nrow=1) +
  guides(fill = guide_legend(override.aes = list(size=3))) +
  labs(x="UMAP-1", y="UMAP-2") +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  # ggplot_theme_NoAxes() +
  theme(
    strip.text = element_text(size=rel(1.15)),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position  ="none",
    legend.direction = "vertical",
    legend.title=element_blank()
  )

pdf(file.path(io$outdir,sprintf("%s_expr_Peyser2019_split_cells.pdf",dimred.to.plot)), width=9.5, height=7)
print(p)
dev.off()

#################################################
## Plot, color by gene expression, all samples ##
#################################################

dimred.to.plot <- "UMAP"

genes.to.plot <- c("COL1A1","NANOG","EPCAM","ZEB1","LIN28A","A4GALT","CD75","ST6GAL1","CD177","CD24")
genes.to.plot <- unique(marker_genes.dt$gene)

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "Mki67"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    condition = sce$condition,
    expr = logcounts(sce)[i,],
    x = reducedDim(sce, dimred.to.plot)[,1],
    y = reducedDim(sce, dimred.to.plot)[,2]
  )
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=0.75, stroke=0.01, alpha=0.90, raster.dpi=150) +
    scale_fill_gradient(low = "gray80", high = "red") +
    # facet_wrap(~condition, nrow=2) +
    # scale_fill_viridis() +
    theme_classic() +
    labs(x="UMAP-1", y="UMAP-2", title=i) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust=0.5),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  pdf(file.path(io$outdir,sprintf("%s_expr_Peyser2019_cells.pdf",dimred.to.plot,i)), width=6, height=5)
  print(p)
  dev.off()
}


###################################################
## Plot, color by gene expression, split samples ##
###################################################

dimred.to.plot <- "PCA"

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "Pi16"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce)[i,],
    condition = sce$condition,
    sample = sce$sample,
    x = reducedDim(sce, dimred.to.plot)[,1],
    y = reducedDim(sce, dimred.to.plot)[,2]
  )# %>% merge(cell_metadata.dt,by="cell")
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=1.50, stroke=0.01, alpha=0.90, raster.dpi=150) +
    facet_wrap(~condition, nrow=1) +
    # scale_fill_gradient(low = "gray80", high = "red") +
    scale_fill_viridis() +
    theme_classic() +
    labs(x="UMAP-1", y="UMAP-2", title=i) +
    # ggplot_theme_NoAxes() +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.0), color="black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust=0.5),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  pdf(file.path(io$outdir,sprintf("%s_expr_cell_cycle_Peyser2019_fibroblasts.pdf",dimred.to.plot,i)), width=9.5, height=7)
  print(p)
  dev.off()
}


##########
## Milo ##
##########

suppressPackageStartupMessages(library(miloR))

# Create Milo object
milo.obj <- Milo(sce_filt[,sce_filt$condition_rename%in%c("Saline_Control","Bleomycin_day11")])

# Build a k-nearest neighbour graph
milo.obj <- buildGraph(milo.obj, k = 20, d = 5, reduced.dim = "PCA")

## Define neighbourhoods on a graph
milo.obj <- makeNhoods(milo.obj, prop = 0.5, k = 15, d=5, refined = TRUE, reduced_dims = "PCA")

# Plot neighbourhood size
plotNhoodSizeHist(milo.obj)

# Counting cells in neighbourhoods
milo.obj <- countCells(milo.obj, meta.data = as.data.frame(colData(milo.obj)), sample="sample")
head(nhoodCounts(milo.obj))

# Calculate within neighbourhood distances
milo.obj <- calcNhoodDistance(milo.obj, d=3, reduced.dim = "PCA", use.assay = "logcounts")

# Defining experimental design
design_matrix.df <- data.frame(colData(milo.obj)[,c("sample","condition_rename")]) %>% unique
design_matrix.df$condition <- as.character(design_matrix.df$condition)
colnames(design_matrix.df) <- c("sample","condition")
rownames(design_matrix.df) <- design_matrix.df$sample
design_matrix.df <- design_matrix.df[colnames(milo.obj@nhoodCounts),,drop=F]

stopifnot(rownames(design_matrix.df)==colnames(milo.obj@nhoodCounts))

DA_results <- testNhoods(milo.obj, design = ~condition, design.df = design_matrix.df[,c("condition"),drop=F], reduced.dim = "PCA")
# contrasts <- "Bleomycin_day11 - Saline_Control"
# DA_results <- testNhoods(milo.obj, design = ~condition, design.df = design_matrix.df[,c("condition"),drop=F], model.contrasts = contrasts, reduced.dim = "PCA")

# Overlay Milo results on the dimensionality reduction
milo.obj <- buildNhoodGraph(milo.obj)

# Plot
p1 <- plotReducedDim(milo.obj, dimred = "PCA", colour_by="condition", text_by = "condition", text_size = 3, point_size=1.50) +
  guides(fill="none")
p2 <- plotNhoodGraphDA(milo.obj, DA_results, layout="PCA",alpha=0.1) 
p1 + p2 + patchwork::plot_layout(guides="collect")

################
## Clustering ##
################

library(bluster)
clusters.test <- clusterCells(sce_filt, use.dimred="PCA", BLUSPARAM=SNNGraphParam(k=50, type="rank", cluster.fun="walktrap"))
table(clusters.test)

sce_filt$cluster <- clusters.test

plotReducedDim(sce_filt, dimred="PCA", colour_by="cluster", point_size=1, point_alpha=1)


#############################
## Save dimred coordinates ##
#############################

# pca.dt <- reducedDim(sce, "PCA") %>% round(3) %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")
# fwrite(pca.dt, file.path(io$outdir,"precomputed_pca.txt.gz"), sep="\t", quote=F)
# 
# umap.dt <- reducedDim(sce, "UMAP") %>% round(3)  %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")
# fwrite(umap.dt, file.path(io$outdir,"precomputed_umap.txt.gz"), sep="\t", quote=F)
# 
# pca_per_sample.dt <- reducedDim(sce, "PCA_per_sample") %>% round(3)  %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")
# fwrite(pca_per_sample.dt, file.path(io$outdir,"precomputed_pca_per_sample.txt.gz"), sep="\t", quote=F)
# 
# umap_per_sample.dt <- reducedDim(sce, "UMAP_per_sample") %>% round(3)  %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")
# fwrite(umap_per_sample.dt, file.path(io$outdir,"precomputed_umap_per_sample.txt.gz"), sep="\t", quote=F)
