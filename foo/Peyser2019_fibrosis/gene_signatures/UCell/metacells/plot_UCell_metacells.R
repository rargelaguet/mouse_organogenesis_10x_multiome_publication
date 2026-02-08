source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(UCell))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"results/mouse/metacells/per_tissue/metacells_metadata.txt.gz")
io$sce <- file.path(io$basedir,"results/mouse/metacells/per_tissue/SingleCellExperiment_metacells.rds")
io$outdir <- paste0(io$basedir,"/results/mouse/gene_signatures/metacells"); dir.create(io$outdir, showWarnings = F)
io$outdir.boxplots <- paste0(io$basedir,"/results/mouse/gene_signatures/metacells/boxplots"); dir.create(io$outdir.boxplots, showWarnings = F)
io$outdir.dimred <- paste0(io$basedir,"/results/mouse/gene_signatures/metacells/dimred"); dir.create(io$outdir.dimred, showWarnings = F)

# Options
opts$tissues <- c(
  "Bone", "Visceral Adipose", "Omentum", "Intestine", "Heart", 
  "Subcutaneous Adipose", "Liver", "Lymph Node", "Lung", "Mesentery", 
  "Muscle", "Pancreas", "Skin", "Spleen", "Tendon", "Artery"
)

opts$clusters <- c(
  "Ccl19", "Cxcl12", "Comp", "Col15a1", "Coch", "Pi16", "Fbln1", "Npnt", "Hhip", "Bmp4"
)

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[Tissue%in%opts$tissues] %>%
  .[ClustName%in%opts$clusters] %>%
  setnames("metacell","cell")

table(cell_metadata.dt$ClustName)
table(cell_metadata.dt$Tissue)

########################
## Load gene metadata ##
########################

# gene_metadata.dt <- fread(io$gene_metadata)

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = F)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

# Load precomputed UCell scores
# altExp(sce,'UCell') <- readRDS(file.path(io$outdir,"UCell_SummarizedExperiment.rds"))[,cell_metadata.dt$cell]
altExp(sce,'UCell') <- readRDS(file.path(io$outdir,"UCell_SummarizedExperiment_GRN.rds"))[,cell_metadata.dt$cell]

##########################
## Load gene signatures ##
##########################

# source(here::here("gene_signatures/load_gene_signatures.R"))

# grep("reik",names(gene.sets),value=T)
# gene.sets["hdf_reik#chronological_age-Down"]

###################################
## Dimensionality reduction plot ##
###################################

dimred.to.plot <- "UMAP"

gene.sets.to.plot <- rownames(altExp(sce, 'UCell'))
gene.sets.to.plot <- grep("K_regulon",rownames(assay(altExp(sce, 'UCell'))), value=T)
# gene.sets.to.plot <- grep("BJ",rownames(assay(altExp(sce, 'UCell'))), value=T)

# i <- gene.sets.to.plot[1]
for (i in gene.sets.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    Tissue = sce$Tissue,
    ClustName = sce$ClustName,
    signature = assay(altExp(sce, 'UCell'))[i,],
    x = reducedDim(sce, dimred.to.plot)[,1],
    y = reducedDim(sce, dimred.to.plot)[,2]
  )
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=signature), shape=21, size=2, stroke=0.01, alpha=0.90, raster.dpi=150) +
    # facet_wrap(~sample, ncol=6) +
    # scale_fill_gradient(low = "gray80", high = "purple") +
    scale_fill_viridis(option="A") +
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
  
  pdf(file.path(io$outdir.dimred,sprintf("%s_signature_%s_Buechler2021_metacells.pdf",i,dimred.to.plot)), width=7, height=6)
  print(p)
  dev.off()
}

###############
## Box plots ##
###############

# gene.sets.to.plot <- rownames(assay(altExp(sce, 'UCell')))
gene.sets.to.plot <- grep("reik",rownames(assay(altExp(sce, 'UCell'))), value=T)

# i <- "KEGG#Oxidative_Phosphorylation_UCell"
for (i in gene.sets.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    signature = assay(altExp(sce, 'UCell'))[i,]
  ) %>% merge(cell_metadata.dt,by="cell")
  
  # to.plot <- to.plot[,N:=.N,by=c("stage","celltype_liu2020")] %>% .[N>=25]
  
  p <- ggplot(to.plot, aes(x=ClustName, y=signature)) +
    # geom_violin(fill="gray70", color="black") +
    geom_boxplot(outlier.shape=NA, fill="gray80", color="black") +
    # stat_summary(fun.data = give.n, fun.args = list(ymax = max(to.plot$expr)+1), size=2.5, geom = "text", position = position_dodge(width=0.75)) +
    # ggrastr::geom_point_rast(size=0.25, alpha=0.80, raster.dpi=150) +
    # scale_fill_gradient(low = "gray80", high = "purple") +
    # geom_hline(yintercept=0, linetype="dashed") +
    # coord_cartesian(ylim=c(0,max(0.15,max(to.plot$signature)))) +
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
  
  pdf(file.path(io$outdir.boxplots,sprintf("%s_signature_boxplots_Buechler2021_metacells.pdf",i)), width=8, height=6)
  print(p)
  dev.off()
}

