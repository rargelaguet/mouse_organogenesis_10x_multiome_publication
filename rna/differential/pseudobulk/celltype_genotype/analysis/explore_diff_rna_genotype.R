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
  "Rostral_neurectoderm",
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
io$marker_genes <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_genes_filtered.txt.gz")
io$diff_results <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype_genotype/parsed/diff_expr_results.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype_genotype/parsed/pdf"); dir.create(io$outdir, showWarnings = F)

###############
## Load data ##
###############

diff.dt <- fread(io$diff_results) %>% 
  .[celltype%in%opts$celltypes] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>%
  .[,logFC_v2:=mean_groupB-mean_groupA] %>% # temporary
  .[, sign := ifelse(logFC>0,"Upregulated in KO","Downregulated in KO")] %>%
  .[, sig := (padj_fdr<=0.01 & abs(logFC)>=2)]

# Print stats
print(sprintf("Number of celltypes: %s",length(unique(diff.dt$celltype))))
print(sprintf("Number of genes: %s",length(unique(diff.dt$gene))))


####################
## Filter results ##
####################

# Filter out genes manually
# dt.filt <- dt.filt[!grep("[^Rik|^Gm|^Rpl]",gene)]
diff.dt <- diff.dt[!grep("^Hb",gene)]

# Subset to lineage markers
marker_genes.dt <- fread(io$marker_genes)
diff_markers.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

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

# to.plot <- diff.dt[sig==T] %>% .[,.N, by=c("celltype","sign")]
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
  .[,direction:=c("Downregulated in KO","Upregulated in KO")[as.numeric(logFC>0)+1]] %>%
  .[,.N, by=c("celltype","direction")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = direction), color="black", stat="identity") + 
  # facet_wrap(~class, scales="free_y") +
  labs(x="", y="Number of DE genes") +
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

to.plot <- diff_markers.dt %>%
  merge(
    marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE
  ) %>% .[,sum(sig), by=c("celltype","celltype_marker","sign")]

p <- ggplot(to.plot, aes(x=celltype, y=V1)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  facet_wrap(~sign, scales="fixed", nrow=2) +
  # facet_wrap(~class, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype_marker)], drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.75))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_fate_bias.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()

################################
## Cell fate bias polar plots ##
################################

to.plot <- diff_markers.dt %>%
  merge(
    marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE
  ) %>% .[,sum(sig,na.rm=T), by=c("celltype","celltype_marker","sign")]

to.plot[V1>=7,V1:=7]

i <- "NMP"
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

##########
## TEST ##
##########

sce <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/test/results/rna/pseudobulk/celltype_genotype/SingleCellExperiment_pseudobulk_with_replicates.rds")

# parse annotations
if (!"celltype" %in% colnames(colData(sce))) {
  sce$celltype <- colnames(sce) %>% strsplit("-") %>% map_chr(1)
}
if (!"genotype" %in% colnames(colData(sce))) {
  sce$genotype <- colnames(sce) %>% strsplit("-") %>% map_chr(2) %>% strsplit("_rep") %>% map_chr(1)
}
if (!"celltype_genotype" %in% colnames(colData(sce))) {
  sce$celltype_genotype <- sprintf("%s_%s",sce$celltype,sce$genotype)
}

tmp <- data.table(sample=colnames(sce), celltype=sce$celltype, genotype=sce$genotype)

genes.to.plot <- c("Fgf8")

to.plot <- as.matrix(logcounts(sce[genes.to.plot,])) %>% 
  as.data.table(keep.rownames = "gene") %>% 
  melt(id.vars = "gene", variable.name = "sample", value.name = "expr") %>% 
  merge(tmp, by="sample")

to.plot[,expr:=expr+0.1]

ggplot(to.plot, aes(x=celltype, y=expr)) +
  geom_bar(aes(fill=genotype), stat="identity", position="dodge", data=to.plot[,.(expr=mean(expr)),by=c("celltype","genotype","gene")]) +
  geom_point(aes(fill=genotype), shape=21, color="black") +
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
