source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(UCell))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(BiocParallel))

set.seed(42)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/rna/scTE/SingleCellExperiment.rds")
io$outdir <- paste0(io$basedir,"/results/rna/gene_signatures/cells"); dir.create(io$outdir, showWarnings = F, recursive = T)
io$outdir.boxplots <- paste0(io$basedir,"/results/rna/gene_signatures/cells/boxplots"); dir.create(io$outdir.boxplots, showWarnings = F)
io$outdir.dimred <- paste0(io$basedir,"/results/rna/gene_signatures/cells/dimred"); dir.create(io$outdir.dimred, showWarnings = F)

# Options
opts$samples <- c( 
  "2F_8OFF",
  "2F_6OFF-2ON",
  "2F_4OFF-4ON",
  "2F_2OFF-2ON-4OFF",
  "2F_4ON-4OFF",
  "2F_8ON",
  "4F_8OFF",
  "4F_6OFF-2ON",
  "4F_4OFF-4ON",
  "4F_2OFF-2ON-4OFF",
  "4F_4ON-4OFF",
  "4F_8ON"
)
opts$cores <- 8

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[,sample:=NULL] %>% setnames("alias","sample") %>% 
  .[sample%in%opts$samples] %>% .[,sample:=factor(sample,levels=opts$samples)]

table(cell_metadata.dt$sample)

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

##########################
## Load gene signatures ##
##########################

source(here::here("rna/plot_gene_signatures/load_gene_signatures.R"))


###############
## Load data ##
###############

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = T)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

# Load precomputed UCell scores
altExp(sce,'UCell') <- readRDS(file.path(io$outdir,"UCell_SummarizedExperiment.rds"))[,colnames(sce)]

# Load pre-computed dimred
io$dimred <- file.path(io$basedir,"results/rna/dimensionality_reduction/cells/2F")
umap.mtx <- fread(file.path(io$dimred,"precomputed_umap.txt.gz")) %>% matrix.please
stopifnot(colnames(sce)%in%rownames(umap.mtx))
umap.mtx <- umap.mtx[colnames(sce),]
reducedDim(sce, "UMAP") <- umap.mtx


###################################
## Dimensionality reduction plot ##
###################################

dimred.to.plot <- "UMAP"

gene.sets.to.plot <- rownames(assay(altExp(sce, 'UCell')))
gene.sets.to.plot <- grep("BJ",rownames(assay(altExp(sce, 'UCell'))), value=T)
gene.sets.to.plot <- grep("KEGG",rownames(assay(sce, 'UCell')), value=T)

# i <- gene.sets.to.plot[1]
for (i in gene.sets.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    sample = sce$sample,
    cell_line = sce$cell_line,
    signature = assay(altExp(sce, 'UCell'))[i,],
    x = reducedDim(sce, dimred.to.plot)[,1],
    y = reducedDim(sce, dimred.to.plot)[,2]
  )
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    geom_point_rast(aes(fill=signature), shape=21, size=0.74, stroke=0.01, alpha=0.90, raster.dpi=150) +
    facet_wrap(~sample, ncol=3) +
    scale_fill_gradient(low = "gray80", high = "purple") +
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
  
  pdf(file.path(io$outdir.dimred,sprintf("%s_signature_%s_BJ_multiome_cells.pdf",i,dimred.to.plot)), width=8, height=5)
  print(p)
  dev.off()
}

###############
## Box plots ##
###############

gene.sets.to.plot <- rownames(assay(altExp(sce, 'UCell')))
# gene.sets.to.plot <- grep("KEGG",rownames(assay(altExp(sce, 'UCell'))), value=T)

# i <- gene.sets.to.plot[1]
for (i in gene.sets.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    signature = assay(altExp(sce, 'UCell'))[i,]
  ) %>% merge(cell_metadata.dt,by="cell")
  
  p <- ggplot(to.plot, aes(x=sample, y=signature, fill=sample)) +
    geom_violin() +
    geom_boxplot(outlier.shape=NA) +
    # stat_summary(fun.data = give.n, fun.args = list(ymax = max(to.plot$expr)+1), size=2.5, geom = "text", position = position_dodge(width=0.75)) +
    # geom_point_rast(size=0.25, alpha=0.80, raster.dpi=150) +
    # scale_fill_gradient(low = "gray80", high = "purple") +
    # geom_hline(yintercept=0, linetype="dashed") +
    facet_wrap(~cell_line, scales="free_x", nrow=1) +
    coord_cartesian(ylim=c(0,max(0.15,max(to.plot$signature)))) +
    guides(x = guide_axis(angle = 90)) +
    labs(x="", y="UCell score", title=i) +
    theme_bw() +
    # ggplot_theme_NoAxes() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(color="black", size=rel(1)),
      axis.text.x = element_text(color="black", size=rel(0.90)),
      axis.title = element_text(color="black", size=rel(1.10)),
      legend.title = element_blank(),
      legend.position = "none"
    )
  
  pdf(file.path(io$outdir.boxplots,sprintf("%s_signature_boxplots_BJ_multiome_cells.pdf",i)), width=8, height=4)
  print(p)
  dev.off()
}
