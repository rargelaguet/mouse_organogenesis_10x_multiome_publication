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
io$sce <- file.path(io$basedir,"processed/cellranger/SingleCellExperiment.rds")
io$outdir.dimred <- file.path(io$basedir,"results/dimensionality_reduction"); dir.create(io$outdir, showWarnings = F)

# Options
opts$samples <- c( 
  # "hTert_t2",
  # "hTert_t4",
  # "hTert_t5",
  # "hTert_t6",
  # "hTert_t7",
  "PDL_25_rep1",
  "PDL_29_rep1",
  "PDL_29_rep2",
  "PDL_33_rep1",
  "PDL_33_rep2",
  "PDL_37_rep1",
  "PDL_37_rep2",
  "PDL_46_rep1",
  "PDL_46_rep2",
  "PDL_50_rep1",
  "PDL_50_rep2"
)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[,sample:=NULL] %>% setnames("alias","sample") %>%
  .[sample%in%opts$samples] %>% .[,sample:=factor(sample,levels=opts$samples)]# %>%
  # .[,stage:=factor(stage, levels=opts$stages, ordered = T)]

cell_metadata.dt[grep("PDL_50",sample),c("cell_cycle_s_score","cell_cycle_g2m_score","cell_cycle_phase"):=list(NA,NA,"G1")]

table(cell_metadata.dt$sample)
table(cell_metadata.dt$genetic_background)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=3000) %>% rownames
sce_filt <- sce[hvgs,]

#########
## PCA ##
#########

sce_filt <- runPCA(sce_filt, ncomponents=15, ntop=3000)
reducedDims(sce) <- reducedDims(sce_filt)

# Plot PCA
# plotPCA(sce, colour_by="sample", point_size=1, point_alpha=1, ncomponents = c(1,2))
plotPCA(sce, colour_by="TOP2A", point_size=1, point_alpha=1)

##########
## UMAP ##
##########

# sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 25, min_dist = 0.3)
# reducedDims(sce) <- reducedDims(sce_filt)

# rotate UMAP
# reducedDim(sce,"UMAP")[,1] <- -reducedDim(sce,"UMAP")[,1]

plotReducedDim(sce, dimred="UMAP", colour_by="group", point_size=1, point_alpha=1)
plotReducedDim(sce, dimred="UMAP", colour_by="cell_cycle_phase", point_size=1, point_alpha=1)
plotReducedDim(sce, dimred="UMAP", colour_by="LMNA", point_size=1, point_alpha=1)
plotReducedDim(sce, dimred="UMAP", colour_by="POU5F1", point_size=1, point_alpha=1)


############################################################
## Load precomputed dimensionality reduction (per donor) ##
############################################################

# io$dimred.dir <- paste0(io$basedir,"/results/dimensionality_reduction/per_sample")
# 
# pca_per_sample.mtx <- opts$samples %>% map(function(i) fread(sprintf("%s/%s/X_pca.txt",io$dimred.dir,i))) %>% rbindlist %>% matrix.please %>% .[colnames(sce),]
# colnames(pca_per_sample.mtx) <- sprintf("PC_%s",1:ncol(pca_per_sample.mtx))
# 
# umap_per_sample.mtx <- opts$samples %>% map(function(i) fread(sprintf("%s/%s/X_umap.txt",io$dimred.dir,i))) %>% rbindlist %>% matrix.please %>% .[colnames(sce),]
# colnames(umap_per_sample.mtx) <- sprintf("UMAP_%s",1:ncol(umap_per_sample.mtx))
# 
# reducedDim(sce, "PCA_per_sample") <- pca_per_sample.mtx
# reducedDim(sce, "UMAP_per_sample") <- umap_per_sample.mtx

#############################################################
## Load precomputed dimensionality reduction (all samples) ##
#############################################################

# io$dimred.dir <- paste0(io$basedir,"/results/dimensionality_reduction")
# 
# pca.mtx <- fread(sprintf("%s/X_pca.txt.gz",io$dimred.dir)) %>% matrix.please %>% .[colnames(sce),]
# colnames(pca.mtx) <- sprintf("PC_%s",1:ncol(pca.mtx))
# 
# umap.mtx <- fread(sprintf("%s/X_umap.txt.gz",io$dimred.dir)) %>% matrix.please %>% .[colnames(sce),]
# colnames(umap.mtx) <- sprintf("UMAP_%s",1:ncol(umap.mtx))
# 
# reducedDim(sce, "PCA") <- pca.mtx
# reducedDim(sce, "UMAP") <- umap.mtx

###################################
## Plot, colour by cell metadata ##
###################################

# dimred.to.plot <- "UMAP_per_sample"
# dimred.to.plot <- "PCA_per_sample"
dimred.to.plot <- "PCA"
# dimred.to.plot <- "UMAP"

dims.to.plot <- c("PC1","PC2")

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

# to.plot[condition=="base",condition:="-FGF2"]
# to.plot[condition=="rich",condition:="+FGF2"]

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  # ggrastr::geom_point_rast(aes(fill=cell_line), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  # ggrastr::geom_point_rast(aes(color=condition), size=1.20, stroke=0.05, alpha=0.90, raster.dpi=300) +
  ggrastr::geom_point_rast(aes(color=cell_cycle_phase), size=0.75, alpha=0.90, raster.dpi=300) +
  guides(color = guide_legend(override.aes = list(size=2.75)), shape = guide_legend(override.aes = list(size=3))) +
  # ggrepel::geom_text_repel(aes(label=condition), size=4.5, data=to.plot[,.(V1=mean(V1), V2=mean(V2)), by=c("condition")]) +
  ggrepel::geom_text_repel(aes(label=cell_cycle_phase), data=to.plot[,.(V1=mean(V1), V2=mean(V2)), by=c("cell_cycle_phase","condition")]) +
  # facet_wrap(~cell_line, nrow=2) +
  # labs(x="UMAP-1", y="UMAP-2") +
  labs(x=sprintf("%s (%.2f%%)",dims.to.plot[1], var.explained[dims.to.plot[1]]), y=sprintf("%s (%.2f%%)",dims.to.plot[2], var.explained[dims.to.plot[2]])) +
  # scale_shape_manual(values=c(21,22)) +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  # ggplot_theme_NoAxes() +
  theme(
    strip.text = element_text(size=rel(1.10)),
    legend.text = element_text(size=rel(1.10)),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    legend.position  ="none",
    # legend.direction = "vertical",
    legend.title=element_blank()
  )

pdf(file.path(io$outdir,sprintf("%s_sample_HDF_media_reik_lab_cells.pdf",dimred.to.plot)), width=12, height=7)
print(p)
dev.off()


##################################################
## Plot, colour by cell metadata, split samples ##
##################################################

dimred.to.plot <- "PCA"

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  # ggrastr::geom_point_rast(aes(fill=condition), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  ggrastr::geom_point_rast(aes(fill=cell_cycle_phase), shape=21, size=1.30, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  guides(fill = guide_legend(override.aes = list(size=3))) +
  # ggrepel::geom_text_repel(aes(label=cell_cycle_phase), data=to.plot[,.(V1=mean(V1), V2=mean(V2)), by=c("group","cell_cycle_phase")]) +
  facet_wrap(~group, nrow=2) +
  labs(x="PC1", y="PC2") +
  theme_classic() +
  theme(
    strip.text = element_text(size=rel(1.10)),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    legend.position  ="right",
    legend.direction = "vertical",
    legend.title=element_blank()
  )

# pdf(file.path(io$outdir,sprintf("%s_cells_split_HDF_media_reik.pdf",dimred.to.plot)), width=10, height=7)
pdf(file.path(io$outdir,sprintf("%s_cells_split_HDF_media_reik.pdf",dimred.to.plot)), width=10, height=8)
print(p)
dev.off()


####################################
## Plot, color by gene expression ##
####################################

dimred.to.plot <- "UMAP"

genes.to.plot <- c("COL1A1","NANOG","EPCAM","ZEB1","LIN28A","A4GALT","CD75","ST6GAL1","CD177","CD24")
# genes.to.plot <- c("DUXA", "ZSCAN4", "MBD3L2","TRIM43", "RFPL1","LEUTX")
genes.to.plot <- marker_genes.dt[group=="pr4"& score>=0.85,gene] %>% head(n=25)
genes.to.plot <- unique(marker_genes.dt$gene)
genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "IL6"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = assay(sce,"logcounts")[i,],
    cell_line = sce$cell_line,
    condition = sce$condition,
    x = reducedDim(sce, dimred.to.plot)[,1],
    y = reducedDim(sce, dimred.to.plot)[,2]
  )
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=0.75, stroke=0.01, alpha=0.90, raster.dpi=150) +
    scale_fill_gradient(low = "gray80", high = "red") +
    # scale_fill_viridis() +
    # facet_wrap(~condition, nrow=2) +
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
  
  pdf(file.path(io$outdir.expr,sprintf("%s_expr_%s_HDF_media_reik_lab_cells.pdf",i,dimred.to.plot)), width=7, height=6)
  print(p)
  dev.off()
}


######################################################
## Plot, color by gene expression, split per sample ##
######################################################

dimred.to.plot <- "PCA_per_sample"

genes.to.plot <- c("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN","BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CD44","CD59","CDH11","CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2","COMP","COPA","CRLF1","CTGF","CTHRC1","CXCL1","CXCL12","CXCL6","CYR61","DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN","EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2","FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A","GADD45B","GAS1","GEM","GJA1","GLIPR1","GLT25D1","GPC1","GPX7","GREM1","HTRA1","ID2","IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","IL8","INHBA","ITGA2","ITGA5","ITGAV","ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","LEPRE1","LGALS1","LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7","MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK","NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","PCOLCE","PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1","SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG","SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2","TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC","TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A")
genes.to.plot <- c("POU5F1", "SOX2", "MYC", "NANOG", "PODXL", "SALL4", "FUT4", "FUT3", "UTF1", "LY6D", "CD9", "NANOS1", "MMD", "CA2", "EPCAM", "ECOR2","DUSP4", "FGF5", "ADCY2", "BRINP1", "DPPA2", "DPPA3", "DPPA5", "KLF4", "DNMT3L", "PRDM14", "TFCP2L1", "KLF17")
genes.to.plot <- marker_genes_liu2020.dt[group=="fm1"& score>=0.85,gene] %>% head(n=50)
genes.to.plot <- unique(marker_genes.dt$gene)

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

# i <- "COL1A1"
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce)[i,],
    cell_line = sce$cell_line,
    condition = sce$condition,
    sample = sce$sample,
    x = reducedDim(sce, dimred.to.plot)[,1],
    y = reducedDim(sce, dimred.to.plot)[,2]
  )# %>% merge(cell_metadata.dt,by="cell")
  
  p <- ggplot(to.plot, aes(x=x, y=y)) +
    ggrastr::geom_point_rast(aes(fill=expr), shape=21, size=0.75, stroke=0.01, alpha=0.90, raster.dpi=150) +
    facet_wrap(~sample, nrow=2) +
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
  
  pdf(file.path(io$outdir.expr,sprintf("%s_expr_%s_HDF_media_reik_lab_63F_cells.pdf",i,dimred.to.plot)), width=12, height=7)
  print(p)
  dev.off()
}


#########################################
## Save precomputed dimred coordinates ##
#########################################

# pca.dt <- reducedDim(sce, "PCA") %>% round(3) %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")
# fwrite(pca.dt, file.path(io$outdir,"precomputed_pca.txt.gz"), sep="\t", quote=F)
# 
# umap.dt <- reducedDim(sce, "UMAP") %>% round(3)  %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")
# fwrite(umap.dt, file.path(io$outdir,"precomputed_umap.txt.gz"), sep="\t", quote=F)

# umap_3d.dt <- reducedDim(sce, "UMAP_3D") %>% round(3)  %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")
# fwrite(umap_3d.dt, file.path(io$outdir,"precomputed_umap_3D.txt.gz"), sep="\t", quote=F)

