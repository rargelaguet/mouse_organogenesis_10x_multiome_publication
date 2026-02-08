# here::i_am("dimensionality_reduction/dimensionality_reduction_sce.R")

library(scran)
library(scater)

source(here::here("settings.R"))
source(here::here("utils.R"))

set.seed(42)

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction"); dir.create(io$outdir, showWarnings = F)
io$outdir.expr <- paste0(io$basedir,"/results/gene_expression/cells"); dir.create(io$outdir.expr, showWarnings = F)

# Options
opts$conditions <- c("Saline Control", "Bleomycin Fibrosis Model", "Bleomycin Fibrosis plus Nintedanib")

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] %>% .[,condition:=factor(condition,levels=opts$conditions)]

table(cell_metadata.dt$condition)

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

dimred.to.plot <- "X_umap"

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  # ggrastr::geom_point_rast(aes(fill=cell_cycle_phase), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  ggrastr::geom_point_rast(aes(fill=cell_type), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  ggrepel::geom_text_repel(aes(label=cell_type), data=to.plot[,.(V1=median(V1), V2=median(V2)), by=c("cell_type")]) +
  guides(fill = guide_legend(override.aes = list(size=3))) +
  labs(x="UMAP-1", y="UMAP-2") +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  # ggplot_theme_NoAxes() +
  theme(
    strip.text = element_text(size=rel(1.15)),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    # legend.position  ="none",
    legend.position  ="right",
    legend.direction = "vertical",
    legend.title=element_blank()
  )

pdf(file.path(io$outdir,sprintf("2F_%s_cell_cycle_Peyser2019_cells.pdf",dimred.to.plot)), width=6, height=5)
print(p)
dev.off()

############################
## Plot, split conditions ##
############################

dimred.to.plot <- "X_umap"

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(fill=cell_type), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  ggrepel::geom_text_repel(aes(label=cell_type), size=2.5,data=to.plot[,.(V1=median(V1), V2=median(V2)), by=c("cell_type","condition")]) +
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

pdf(file.path(io$outdir,sprintf("%s_celltype_Peyser2019_split_cells.pdf",dimred.to.plot)), width=10, height=6)
print(p)
dev.off()

#########################
## Plot, split samples ##
#########################

dimred.to.plot <- "X_umap"

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(fill=cell_type), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  ggrepel::geom_text_repel(aes(label=cell_type), size=2.5, data=to.plot[,.(V1=median(V1), V2=median(V2)), by=c("cell_type","sample")]) +
  # facet_wrap(~condition, nrow=1) +
  facet_wrap(~sample, nrow=3) +
  guides(fill = guide_legend(override.aes = list(size=3))) +
  labs(x="UMAP-1", y="UMAP-2") +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  # ggplot_theme_NoAxes() +
  theme(
    strip.text = element_text(size=rel(0.65)),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position  ="none",
    legend.direction = "vertical",
    legend.title=element_blank()
  )

pdf(file.path(io$outdir,sprintf("%s_celltype_Peyser2019_split_cells.pdf",dimred.to.plot)), width=10, height=6)
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

dimred.to.plot <- "UMAP"

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "MKI67"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce)[i,],
    condition = sce$condition,
    cell_line = sce$cell_line,
    sample = sce$sample,
    x = reducedDim(sce, dimred.to.plot)[,1],
    y = reducedDim(sce, dimred.to.plot)[,2]
  )# %>% merge(cell_metadata.dt,by="cell")
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=0.75, stroke=0.01, alpha=0.90, raster.dpi=150) +
    # facet_wrap(~condition, nrow=2) +
    facet_wrap(~cell_line, nrow=1) +
    scale_fill_gradient(low = "gray80", high = "red") +
    # scale_fill_viridis() +
    theme_classic() +
    labs(x="UMAP-1", y="UMAP-2", title=i) +
    # ggplot_theme_NoAxes() +
    theme_classic() +
    theme(
      strip.text = element_text(size=rel(1.15), color="black"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust=0.5),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  pdf(file.path(io$outdir,sprintf("2F_%s_expr_cell_cycle_Peyser2019_cells.pdf",dimred.to.plot,i)), width=9.5, height=7)
  print(p)
  dev.off()
}

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
