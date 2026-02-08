# here::i_am("dimensionality_reduction/dimensionality_reduction.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
io$sce <- file.path(io$basedir,"processed/atac/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/atac/dimensionality_reduction"); dir.create(io$outdir, showWarnings = F)

# Options
opts$groups <- c(
  "p8_Control", 
  "p11_Control", 
  "p14_Control", 
  "p19_Control", 
  "p22_Control", 
  "p24_Control", 
  "p29_Control", 
  "p8_TG", 
  "p11_TG", 
  "p14_TG", 
  "p19_TG", 
  "p22_TG", 
  "p24_TG", 
  "p29_TG"
)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata) %>%
  .[Group%in%opts$groups]

table(sample_metadata.dt$Group)

###############
## Load data ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]

# Add sample metadata as colData
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

# Normalise
assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(assay(sce,"counts")),"/"))+1)

##################
## Subset genes ##
##################

# genes.to.use <- grep("^RP[L|S]",rownames(sce),value=T)
# genes.to.use <- grep("^MT-",rownames(sce),value=T)
# sce.filt <- sce[genes.to.use,]

#######################
## Feature selection ##
#######################

# decomp <- modelGeneVar(sce)
# hvgs <- rownames(decomp)[decomp$p.value<0.01 & decomp$mean > 1]
# sce.filt <- sce[hvgs,]
# dim(sce.filt)

tmp <- apply(assay(sce,"logcounts"),1,var)
sce.filt <- sce[tmp>0.05,]
sce.filt

#########
## PCA ##
#########

sce.filt <- runPCA(sce.filt, ncomponents=5)
reducedDims(sce) <- reducedDims(sce.filt)

# cor(sample_metadata.dt$age,reducedDim(sce.filt,"PCA"))

plotPCA(sce.filt, colour_by="Passage", shape_by="Treatment", ncomponents=c(1,2), point_alpha=1, point_size=5)
plotPCA(sce.filt, colour_by="chr4:10761757-10762940", ncomponents=c(1,2), point_alpha=1, point_size=5)
# plotPCA(sce.filt, colour_by="Passage", ncomponents=c(1,2), point_alpha=1, point_size=2.5)
# plotPCA(sce.filt, colour_by="Group", ncomponents=c(1,2), point_alpha=1, point_size=2.5)

var.explained <- attr(reducedDim(sce.filt,"PCA"),"varExplained")
var.explained <- (100*(var.explained/sum(var.explained))) %>% round(1)

##########
## Plot ##
##########

to.plot <- reducedDim(sce.filt) %>% as.data.table(keep.rownames = T) %>% 
  setnames("rn","sample") %>%
  # melt(id.vars="sample", variable.name="pc") %>%
  merge(sample_metadata.dt,by="sample")

p <- ggplot(to.plot, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill=Passage, shape=Treatment), shape=21, size=5, stroke=0.35) +
  # ggrepel::geom_text_repel(aes_string(label="surface_marker"), data=to.plot[,.(PC1=median(PC1), PC2=median(PC2)), by="surface_marker"]) +
  scale_fill_viridis_c() +
  theme_classic() +
  labs(x=sprintf("PC1 (%.2f%%)",var.explained[1]), y=sprintf("PC1 (%.2f%%)",var.explained[2])) +
  # ggplot_theme_NoAxes() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color="black")
    # axis.ticks = element_blank()
  )

pdf(file.path(io$outdir,"pca_Chan2022_senescence.pdf"), width=8, height=6)
print(p)
dev.off()

######################
## Plot PCA weights ##
######################

pca <- irlba::prcomp_irlba(t(logcounts(sce.filt)))
# pca_gene_weights.dt <- attr(reducedDims(sce.filt)[["PCA"]],"rotation") %>% as.data.table(keep.rownames = T) %>% setnames("rn","gene")
pca_gene_weights.dt <- pca$rotation %>% as.data.table() %>% .[,gene:=rownames(sce.filt)]

# attr(reducedDims(sce)[["PCA"]],"rotation")[,1] %>% abs %>% sort(decreasing = T) %>% head

# plotPCA(sce.filt, colour_by="KRT6A", ncomponents=c(1,2), point_alpha=1, point_size=2.5)

#########################
## GSEA on PCA weights ##
#########################

library(clusterProfiler)
library(org.Hs.eg.db)

symbol2ensembl.dt <- bitr(pca_gene_weights.dt$gene, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db", drop=T) %>%
  as.data.table(keep.rownames=F) %>% setnames(c("gene","ensembl")) %>% .[!(duplicated(gene) | duplicated(ensembl))]


to_plot_pca_gene_weights.dt <- pca_gene_weights.dt  %>% merge(symbol2ensembl.dt, by="gene")
tmp <- to_plot_pca_gene_weights.dt$PC1; names(tmp) <- to_plot_pca_gene_weights.dt$gene
tmp <- sort(tmp,decreasing = T)

# gseaKEGG <- gseKEGG(
#   geneList    = sort(tmp), # ordered named vector of fold changes (Entrez IDs are the associated names)
#   OrgDb       = org.Hs.eg.db,
#   keyType     = 'ENSEMBL',
#   nPerm       = 1000,
#   minGSSize   = 20,
#   pvalueCutoff = 0.01,
#   verbose = FALSE
# )

gseaGO <- gseGO(
  geneList    = tmp, 
  ont         = "MF",
  OrgDb       = org.Hs.eg.db,
  keyType     = 'SYMBOL',
  minGSSize   = 20,
  pvalueCutoff = 0.01
)

View(gseaGO@result)
View(gseaGO@result[,1:8])

  
# library(msigdbr)
# msigdbr_show_species()
# C3_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
#   dplyr::select(gs_name, entrez_gene)

###################################
## PCA using a specific gene set ##
###################################

genes.to.use <- c("HSPA1A", "HSPA1B", "PDLIM1", "IGFBP3", "CLIC1", "IQGAP1", "CXCR2", "TP53", "PLEC", "HMGA1", "CDKN1A", "CDKN2A", "IL6", "ACTB", "CXCL1", "NFKB1", "RELA", "STAT3", "CDKN1A", "CCL2")
sce.small <- runPCA(sce[genes.to.use], ncomponents=2)

var.explained <- attr(reducedDims(sce.small)[["PCA"]],"percentVar") %>% round(2)
# var.explained <- 100*(var.explained / sum(var.explained))

to.plot <- reducedDim(sce.small,"PCA")[,c(1,2)] %>% 
  as.data.table(keep.rownames = T) %>% setnames("rn","sample") %>%
  merge(sample_metadata.dt,by="sample")

p <- ggplot(to.plot, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill=Passage, shape=Treatment), size=5, stroke=0.35) +
  labs(x=sprintf("PC1 (%.2f%%)",var.explained[1]), y=sprintf("PC2 (%.2f%%)",var.explained[2])) +
  scale_fill_viridis_c() +
  scale_shape_manual(values=c(21,22)) +
  # scale_shape_manual(values=c(23,24)) +
  # scale_color_manual(values=brewer.pal(n=9,"OrRd")) +
  theme_classic() +
  theme(
    axis.text = element_text(colour="black",size=rel(0.75)),
    legend.position = "right"
  )

pdf(file.path(io$outdir,"pca_senescence_Chan2022_senescence.pdf"), width=8, height=6)
print(p)
dev.off()

for (i in genes.to.use) {
  
  p <- plotPCA(sce.small, color_by=i, point_size=3, point_alpha=1)
  
  pdf(file.path(io$outdir,sprintf("pca_senescence_Chan2022_senescence_%s.pdf",i)), width=8, height=6)
  print(p)
  dev.off()
}
  
  
  
