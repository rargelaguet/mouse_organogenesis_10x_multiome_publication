source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))

set.seed(42)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"results/metacells/per_tissue/metacells_metadata.txt.gz")
io$sce <- file.path(io$basedir,"results/metacells/per_tissue/SingleCellExperiment_metacells.rds")
io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction/metacells"); dir.create(io$outdir, showWarnings = F)

# Options
opts$tissues <- c(
  "Bone", "Visceral Adipose", "Omentum", "Intestine", "Heart", 
  "Subcutaneous Adipose", "Liver", "Lymph Node", "Lung", "Mesentery", 
  "Muscle", "Pancreas", "Skin", "Spleen", "Tendon", "Artery"
)

opts$clusters <- c(
  "Ccl19", "Cxcl12", "Comp", "Col15a1", "Coch", "Pi16", "Fbln1", "Npnt", "Hhip", "Bmp4"
)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[Tissue%in%opts$tissues] %>%
  .[ClustName%in%opts$clusters] %>%
  setnames("metacell","cell")

table(cell_metadata.dt$ClustName)
table(cell_metadata.dt$Tissue)

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = T)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

table(sce$sample)

####################
## Load gene sets ##
####################

# marker_genes_liu2020.dt <- fread("/Users/rargelaguet/data/Liu2020_human_reprogramming/results/rna/differential/cells/snClust/parsed/marker_genes_filtered.txt.gz") %>%
#   .[gene%in%rownames(sce)]

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=1500) %>% rownames
hvgs <- hvgs[!hvgs%in%c("POU5F1","KLF4","MYC","SOX2")]
sce_filt <- sce[hvgs,]

#########
## PCA ##
#########

sce_filt <- runPCA(sce_filt, ncomponents = 25, ntop=1500)
reducedDims(sce) <- reducedDims(sce_filt)

##########
## UMAP ##
##########

sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 25, min_dist = 0.50)
reducedDims(sce) <- reducedDims(sce_filt)

# plotReducedDim(sce, dimred = "UMAP", colour_by="sample", point_size=2, point_alpha=1)
# plotReducedDim(sce, dimred = "UMAP", colour_by="COL1A1", point_size=2, point_alpha=1)
# plotReducedDim(sce, dimred = "UMAP", colour_by="POU5F1", point_size=2, point_alpha=1)


######################
## Plot all samples ##
######################

dimred.to.plot <- "UMAP"

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(fill=Tissue), shape=21, size=2.40, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  guides(fill = guide_legend(override.aes = list(size=3))) +
  # scale_fill_brewer(palette="Dark2") +
  # scale_fill_manual(values=opts$stage.colors) +
  labs(x="UMAP-1", y="UMAP-2") +
  theme_classic() +
  # ggplot_theme_NoAxes() +
  theme(
    strip.text = element_text(size=rel(1.15)),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position  ="right",
    legend.direction = "vertical",
    legend.title=element_blank()
  )

pdf(file.path(io$outdir,sprintf("%s_tissue_metacells_Buechler2021.pdf",dimred.to.plot)), width=8, height=5)
print(p)
dev.off()

#########################
## Plot, split tissues ##
#########################

dimred.to.plot <- "UMAP"

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(fill=ClustName), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  facet_wrap(~Tissue, nrow=3) +
  guides(fill = guide_legend(override.aes = list(size=3))) +
  # scale_color_viridis() +
  # scale_fill_brewer(palette="Dark2") +
  # ggrepel::geom_text_repel(aes(label=stage), data=to.plot[,.(V1=median(V1), V2=median(V2)), by=c("stage","dataset")]) +
  labs(x="UMAP-1", y="UMAP-2") +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  # ggplot_theme_NoAxes() +
  theme(
    strip.text = element_text(size=rel(1.0)),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position  ="right",
    legend.direction = "vertical",
    legend.title=element_blank()
  )

pdf(file.path(io$outdir,sprintf("%s_tissue_metacells_Buechler2021.pdf",dimred.to.plot)), width=8, height=8)
print(p)
dev.off()


#################################################
## Plot, color by gene expression, all samples ##
#################################################

dimred.to.plot <- "UMAP"

genes.to.plot <- c("COL1A1","NANOG","EPCAM","ZEB1","LIN28A","A4GALT","CD75","ST6GAL1","CD177","CD24")

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "Pi16"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    Tissue = sce$Tissue,
    ClustName = sce$ClustName,
    expr = logcounts(sce)[i,],
    x = reducedDim(sce, dimred.to.plot)[,1],
    y = reducedDim(sce, dimred.to.plot)[,2]
  )
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=1.60, stroke=0.01, alpha=0.90, raster.dpi=150) +
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
  
  pdf(file.path(io$outdir,sprintf("2F_%s_expr_cell_cycle_BJ_multiome_CHAMP_cells.pdf",dimred.to.plot,i)), width=6, height=5)
  print(p)
  dev.off()
}


###################################################
## Plot, color by gene expression, split samples ##
###################################################

dimred.to.plot <- "UMAP"

genes.to.plot <- c("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COMP","COPA","CRLF1","CTGF","CTHRC1","CXCL1","CXCL12","CXCL6","CYR61","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","GLT25D1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","IL8","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","LEPRE1","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A")
genes.to.plot <- c("POU5F1", "SOX2", "MYC", "NANOG", "PODXL", "SALL4", "FUT4", "FUT3", "UTF1", "LY6D", "CD9", "NANOS1", "MMD", "CA2", "EPCAM", "ECOR2","DUSP4", "FGF5", "ADCY2", "BRINP1", "DPPA2", "DPPA3", "DPPA5", "KLF4", "DNMT3L", "PRDM14", "TFCP2L1", "KLF17")
genes.to.plot <- marker_genes_liu2020.dt[group=="fm1"& score>=0.85,gene] %>% head(n=50)
genes.to.plot <- unique(marker_genes.dt$gene)

genes.to.plot <- c("HERVH-int","HERVK-int","HERVL-int","LTR5_Hs","LTR5A","LTR5B","LTR7A","LTR7B","MLT2A1")

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "HAS2"
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
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=2, stroke=0.01, alpha=0.90, raster.dpi=150) +
    facet_wrap(~condition, nrow=2) +
    # scale_fill_gradient(low = "gray80", high = "red") +
    scale_fill_viridis() +
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
  
  pdf(file.path(io$outdir,sprintf("2F_%s_expr_cell_cycle_BJ_multiome_CHAMP_cells.pdf",dimred.to.plot,i)), width=9.5, height=7)
  print(p)
  dev.off()
}
