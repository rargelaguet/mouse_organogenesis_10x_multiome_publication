# here::i_am("rna/differential/wt_vs_ko/analysis/barplots.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Load utils
# source(here::here("differential/analysis/utils.R"))

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  # "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

##############
## Settings ##
##############

# I/O
io$basedir <- file.path(io$basedir,"test")
io$marker_peaks <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype/PeakMatrix/parsed/markers_filt.txt.gz")
io$diff_results <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype_genotype/PeakMatrix/parsed/diff_results.txt.gz")
io$outdir <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype_genotype/PeakMatrix/parsed/pdf"); dir.create(io$outdir, showWarnings = F)

###############
## Load data ##
###############

diff.dt <- fread(io$diff_results) %>% 
  .[celltype%in%opts$celltypes] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>%
  .[,logFC:=mean_groupB-mean_groupA] %>% # temporary
  .[, sign := ifelse(logFC>0,"Upregulated in KO","Downregulated in KO")] %>%
  .[, sig := (padj_fdr<=0.01 & abs(logFC)>=1.5)]

# Print stats
print(sprintf("Number of celltypes: %s",length(unique(diff.dt$celltype))))
print(sprintf("Number of features: %s",length(unique(diff.dt$feature))))

####################
## Filter results ##
####################

# Subset to lineage markers
marker_peaks.dt <- fread(io$marker_peaks)
diff_markers.dt <- diff.dt[feature%in%unique(marker_peaks.dt$feature)]

################
## Polar plot ##
################

to.plot <- diff_markers.dt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype","sign")]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  facet_wrap(~sign, nrow=1) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  # guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
  )

pdf(file.path(io$outdir,"DE_polar_plots_marker_genes.pdf"), width=11, height=8)
print(p)
dev.off()

##############
## Bar plot ##
##############

to.plot <- diff.dt[sig==T] %>% .[,.N, by=c("celltype","sign")]
to.plot <- diff_markers.dt[sig==T] %>% .[,.N, by=c("celltype","sign")]

p <- ggplot(to.plot, aes(x=celltype, y=N, fill=sign)) +
  geom_bar(color="black", stat = 'identity') + 
  # scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "none",
    # axis.line = element_blank(),
    axis.text.y = element_text(color="black", size=rel(1)),
    axis.text.x = element_text(color="black", size=rel(0.9))
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )


pdf(file.path(io$outdir,"DE_barplots_marker_genes.pdf"), width=8, height=4)
print(p)
dev.off()

##########################################
## Barplot of fraction of genes up/down ##
##########################################

to.plot <- diff_markers.dt[sig==T] %>%
  .[,.N, by=c("celltype","sign")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = sign), color="black", stat="identity") + 
  # facet_wrap(~class, scales="free_y") +
  labs(x="", y="Number of DA peaks") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    # axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.75))
  )

pdf(file.path(io$outdir,"DE_barplots_marker_genes_direction.pdf"), width=8, height=4.5)
print(p)
dev.off()

#############################
## Cell fate bias barplots ##
#############################

# to.plot <- diff_markers.dt %>%
#   merge(
#     marker_peaks.dt[,c("feature","celltype")] %>% setnames("celltype","celltype_marker"), by = "feature", allow.cartesian=TRUE
#   ) %>% .[,sum(sig), by=c("celltype","celltype_marker","sign")]
# 
# p <- ggplot(to.plot, aes(x=celltype, y=V1)) +
#   geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
#   facet_wrap(~sign, scales="fixed", nrow=2) +
#   # facet_wrap(~class, scales="fixed") +
#   scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype_marker)], drop=F) +
#   labs(x="", y="Number of DA peaks") +
#   guides(x = guide_axis(angle = 90)) +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     axis.line = element_blank(),
#     axis.text.x = element_text(color="black", size=rel(0.75))
#   )
# 
# pdf(sprintf("%s/DE_barplots_marker_genes_fate_bias.pdf",io$outdir), width=9, height=5)
# print(p)
# dev.off()

################################
## Cell fate bias polar plots ##
################################

to.plot <- diff_markers.dt %>%
  merge(
    marker_peaks.dt[,c("feature","celltype")] %>% setnames("celltype","celltype_marker"), by = "feature", allow.cartesian=TRUE
  ) %>% .[,sum(sig), by=c("celltype","celltype_marker","sign")]

# to.plot[V1>=7,V1:=7]

for (i in unique(to.plot$celltype)) {
  
  p <- ggplot(to.plot[celltype==i], aes(x=celltype_marker, y=V1)) +
    geom_bar(aes(fill = celltype_marker), color="black", stat = 'identity') + 
    facet_wrap(~sign, nrow=1) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    coord_polar() +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.text = element_text(size=rel(0.75)),
      legend.title = element_blank(),
      axis.title=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      axis.text.x = element_blank()
      # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
    )
  
  pdf(file.path(io$outdir,sprintf("%s_DE_polar_plot_fate_bias.pdf",i)), width=6, height=3)
  print(p)
  dev.off()
}

##################
## Volcano plot ##
##################

to.plot <- diff.dt[celltype=="NMP"]
# to.plot <- diff.dt[celltype=="Neural_crest"]

negative_hits <- to.plot[sig==TRUE & logFC<0,feature]
positive_hits <- to.plot[sig==TRUE & logFC>0,feature]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$logFC), na.rm=T)
ylim <- max(-log10(to.plot$padj_fdr), na.rm=T)

to.plot.subset <- rbind(to.plot[sig==TRUE],to.plot[sig==F][sample(x=.N, size=.N/4)])

ggplot(to.plot.subset, aes(x=logFC, y=-log10(padj_fdr))) +
  labs(title="", x="Log Fold Change", y=expression(paste("-log"[10],"(p.value)"))) +
  # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
  geom_segment(x=0, xend=0, y=0, yend=ylim-1, color="orange") +
  ggrastr::geom_point_rast(aes(color=sig), size=1) +
  scale_color_manual(values=c("black","red")) +
  scale_x_continuous(limits=c(-xlim-0.5,xlim+0.5)) +
  scale_y_continuous(limits=c(0,ylim+1)) +
  annotate("text", x=0, y=ylim+1, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim, y=ylim+1, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim, y=ylim+1, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  theme_classic() +
  theme(
    # axis.text=element_text(size=rel(1.75), color='black'),
    # axis.title=element_text(size=rel(1.95), color='black'),
    # axis.title.y = element_text(),
    # axis.title.x = element_text(),
    legend.position="none",
    # panel.border=element_blank(),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank()
    # panel.background = element_blank()
  )

##########
## TEST ##
##########

atac.se <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/test/results/atac/archR/pseudobulk/celltype_genotype/PeakMatrix/PeakMatrix_pseudobulk_with_replicates.rds")

# parse annotations
if (!"celltype" %in% colnames(colData(atac.se))) {
  atac.se$celltype <- colnames(atac.se) %>% strsplit("-") %>% map_chr(1)
}
if (!"genotype" %in% colnames(colData(atac.se))) {
  atac.se$genotype <- colnames(atac.se) %>% strsplit("-") %>% map_chr(2) %>% strsplit("_rep") %>% map_chr(1)
}
if (!"celltype_genotype" %in% colnames(colData(atac.se))) {
  atac.se$celltype_genotype <- sprintf("%s_%s",atac.se$celltype,atac.se$genotype)
}

tmp <- data.table(sample=colnames(atac.se), celltype=atac.se$celltype, genotype=atac.se$genotype)

peaks.to.plot <- c("chr9:56069989-56070589")

to.plot <- as.matrix(logcounts(atac.se[peaks.to.plot,])) %>%
  as.data.table(keep.rownames = "feature") %>%
  melt(id.vars = "feature", variable.name = "sample", value.name = "expr") %>%
  merge(tmp, by="sample")

# to.plot[,expr:=expr+0.1]

ggplot(to.plot, aes(x=celltype, y=expr)) +
  geom_bar(aes(fill=genotype), stat="identity", position="dodge") +
  # facet_wrap(~gene) +
  # scale_fill_gradient(low = "gray80", high = "purple") +
  # labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    # legend.position="none"
  )
