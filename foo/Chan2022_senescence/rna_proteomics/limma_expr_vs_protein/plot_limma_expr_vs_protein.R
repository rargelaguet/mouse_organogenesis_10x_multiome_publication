source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$sample_metadata.expr <- file.path(io$basedir,"processed/rna/sample_metadata.txt")
io$sample_metadata.proteomics <- file.path(io$basedir,"processed/proteomics/sample_metadata.txt")
io$feature_metadata.proteomics <- file.path(io$basedir,"processed/proteomics/feature_metadata.txt")
io$sce.expr <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
io$sce.proteomics <- paste0(io$basedir,"/processed/proteomics/SingleCellExperiment.rds")
io$pca_weights.expr <- paste0(io$basedir,"/results/rna/dimensionality_reduction/PCA_weights.txt.gz")
io$pca_weights.proteomics <- paste0(io$basedir,"/results/proteomics/dimensionality_reduction/PCA_weights.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna_proteomics/rna_vs_protein_levels"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$plot_logcounts_rna <- TRUE

###########################
## Load feature metadata ##
###########################

feature_metadata_proteomics.dt <- fread(io$feature_metadata.proteomics)

######################
## Load PCA weights ##
######################

pca_weights_rna.dt <- fread(io$pca_weights.expr)
pca_weights_proteomics.dt <- fread(io$pca_weights.proteomics) %>% .[,PC1:=-PC1]

################################
## Scatterplot of PCA weights ##
################################

to.plot <- merge(pca_weights_rna.dt[,c("gene","PC2")] %>% setnames("PC2","pca_weight"), pca_weights_proteomics.dt[,c("gene","PC2")] %>% setnames("PC2","pca_weight"), by=c("gene"), suffixes=c("_rna","_proteomics")) %>%
  # .[abs(pca_weight_rna)>0.01 & abs(pca_weight_proteomics)>0.01] %>%
  .[,tmp:=abs(pca_weight_rna)*abs(pca_weight_proteomics)] %>% setorder(-tmp) %>% .[,tmp:=NULL]

ylim.max <- max(abs(to.plot$pca_weight_rna))
xlim.max <- max(abs(to.plot$pca_weight_proteomics))
lim.max <- max(xlim.max,ylim.max)

genes.to.highlight <- c()

p <- ggplot(to.plot, aes(x=pca_weight_rna, y=pca_weight_proteomics)) +
  geom_point(size=0.40, color="gray30") +
  stat_cor(method = "pearson") +
  stat_smooth(method="lm", color="gray70") +
  coord_cartesian(xlim=c(-lim.max,lim.max), ylim=c(-lim.max,lim.max)) +
  geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
  geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
  ggrepel::geom_text_repel(data=to.plot[gene%in%genes.to.highlight], aes(x=pca_weight_rna, y=pca_weight_proteomics, label=gene), color="purple", size=3, segment.size=0.1, max.overlaps = Inf) +
  ggrepel::geom_text_repel(data=to.plot[abs(pca_weight_rna)>=0.025 & abs(pca_weight_proteomics)<=0.025], aes(x=pca_weight_rna, y=pca_weight_proteomics, label=gene), color="purple", size=3, segment.size=0.1, max.overlaps = Inf) +
  ggrepel::geom_text_repel(data=to.plot[abs(pca_weight_proteomics)>=0.025 & abs(pca_weight_rna)<=0.025], aes(x=pca_weight_rna, y=pca_weight_proteomics, label=gene), color="purple", size=3, segment.size=0.1, max.overlaps = Inf) +
  ggrepel::geom_text_repel(data=to.plot[abs(pca_weight_proteomics)>=0.025 & abs(pca_weight_rna)>=0.025], aes(x=pca_weight_rna, y=pca_weight_proteomics, label=gene), color="purple", size=3, segment.size=0.1, max.overlaps = Inf) +
  labs(x="PC2 weights (RNA)", y="PC2 weights (Proteomics)", title="") +
  # facet_wrap(~transcript_type, nrow=1) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(1.0), color="black"),
    legend.position = "none",
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,"scatterplot_PC2_weights_BJ_passaging_multiome_CHAMP.pdf"), width=8, height=7)
print(p)
dev.off()
