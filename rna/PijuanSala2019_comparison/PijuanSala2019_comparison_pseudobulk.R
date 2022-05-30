# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$outdir <- paste0(io$basedir,"/results/rna/PijuanSala2019_comparison"); dir.create(io$outdir, showWarnings = F)
dir.create(file.path(io$outdir,"per_celltype"), showWarnings = F)
dir.create(file.path(io$outdir,"per_gene"), showWarnings = F)

#######################
## Load gene markers ##
#######################

marker_genes.dt <- fread(io$rna.atlas.marker_genes.up)

# Filter out some genes
marker_genes.dt <- marker_genes.dt[grep("^Rik|Rik$|^mt-|^Rps-|^Rpl-|^Gm",gene,invert = T)]

###############################################
## Load pseudobulk expression from the atlas ##
###############################################

# io$average_expression_per_celltype <- paste0(io$basedir,"/results/rna/pseudobulk/avg_expr_per_celltype_and_gene.txt.gz")
sce.atlas.pseuboulk <- readRDS(io$rna.atlas.sce.pseudobulk)[,opts$celltypes]
# colnames(sce.atlas) <- stringr::str_replace_all(colnames(sce.atlas),"[/ ]","_")

# Filter genes
sce.atlas.pseuboulk <- sce.atlas.pseuboulk[rownames(sce.atlas.pseuboulk) %in% unique(marker_genes.dt$ens_id),]

################################
## Load pseudobulk expression ##
################################

# io$average_expression_per_celltype <- paste0(io$basedir,"/results/rna/pseudobulk/avg_expr_per_celltype_and_gene.txt.gz")
sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]
# colnames(sce.atlas) <- stringr::str_replace_all(colnames(sce.atlas),"[/ ]","_")

# Filter genes
sce.pseudobulk <- sce.pseudobulk[rownames(sce.pseudobulk) %in% unique(marker_genes.dt$gene),]

#################################
## Create data.table and merge ##
#################################

pseudobulk_expr.dt <- logcounts(sce.pseudobulk) %>% as.data.table(keep.rownames = T) %>%
  melt(id.vars="rn", variable.name="celltype", value.name="expr") %>% setnames("rn","gene") %>%
  .[,dataset:=as.factor("This study")]

pseudobulk_expr_atlas.dt <- logcounts(sce.atlas.pseuboulk) %>% as.data.table(keep.rownames = T) %>%
  melt(id.vars="rn", variable.name="celltype", value.name="expr") %>% setnames("rn","ens_id") %>%
  merge(unique(marker_genes.dt[,c("ens_id","gene")]),by="ens_id") %>%
  .[,dataset:=as.factor("PijuanSala2019")]

expr.dt <- merge(
  pseudobulk_expr_atlas.dt[,c("celltype","gene","expr")],
  pseudobulk_expr.dt[,c("celltype","gene","expr")],
  by = c("celltype","gene"),
  suffixes = c(".atlas",".query")
)

length(unique(expr.dt$gene))

#####################################
## Calculate correlations per gene ##
#####################################

cor_genes.dt <- expr.dt %>% 
  .[, cor.test(x=expr.atlas, y=expr.query)[c("estimate","p.value")], by = c("gene")] %>%
  # .[, para := rep(c("r","p"), .N/2)] %>% data.table::dcast(gene ~ para, value.var = "V1") %>%
  .[, padj_fdr := p.adjust(p.value, method="fdr")] %>%
  .[, sig := padj_fdr <= 0.1] %>%
  setorder(padj_fdr)

to.plot <- cor_genes.dt

negative_hits <- to.plot[sig==TRUE & estimate<0,gene]
positive_hits <- to.plot[sig==TRUE & estimate>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$estimate), na.rm=T)
ylim <- max(-log10(to.plot$p.value), na.rm=T)

p <- ggplot(to.plot, aes(x=estimate, y=-log10(padj_fdr))) +
  labs(title="", x="Pearson correlation", y=expression(paste("-log"[10],"(q.value)"))) +
  geom_segment(x=0, xend=0, y=0, yend=ylim-1, color="orange") +
  ggrastr::geom_point_rast(aes(color=sig), size=1) +
  scale_color_manual(values=c("black","red")) +
  scale_x_continuous(limits=c(-xlim-0.05,xlim+0.05)) +
  scale_y_continuous(limits=c(0,ylim+1)) +
  annotate("text", x=0, y=ylim+1, size=5, label=sprintf("(%d)", all)) +
  annotate("text", x=-0.75, y=ylim+1, size=5, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=0.75, y=ylim+1, size=5, label=sprintf("%d (+)",length(positive_hits))) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=10), aes(x=estimate, y=-log10(padj_fdr), label=gene), max.overlaps=Inf, size=4) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(1), color='black'),
    axis.title = element_text(size=rel(1), color='black'),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"PijuanSala2019_comparison_volcano_plot.pdf"), width = 6, height = 6)
print(p)
dev.off()

###############################################################
## Explore Cytoplasm/Nuclear expression ratios per cell type ##
###############################################################

# dt[,diff:=expr.atlas-expr.query]
# 
# to.plot <- dt %>% merge(marker_genes.dt[,c("gene","celltype")], by=c("gene","celltype"))
# to.plot <- dt
# 
# order.celltypes <- to.plot[,.(diff=mean(diff)),by="celltype"] %>% setorder(diff) %>% .$celltype
# to.plot[,celltype:=factor(celltype,levels=order.celltypes)]
# 
# p <- ggplot(to.plot, aes(x=celltype, y=diff)) +
#   geom_boxplot(aes(fill=celltype), outlier.shape=NA, alpha=0.8) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   coord_flip() +
#   labs(y="Atlas-Query RNA expression difference") +
#   geom_hline(yintercept=0, linetype="dashed", size=0.5) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     strip.background = element_blank(),
#     strip.text = element_text(color="black", size=rel(1.2)),
#     axis.title.x = element_text(color="black", size=rel(1.1)),
#     axis.title.y = element_blank(),
#     axis.text.y = element_text(size=rel(1), color="black"),
#     axis.text.x = element_text(size=rel(0.75), color="black")
#   )
# 
# pdf(sprintf("%s/nucleus_vs_cell_boxplots.pdf",io$outdir), width = 5, height = 8)
# print(p)
# dev.off()

###########################
## Scatterplots per gene ##
###########################

# genes.to.plot <- unique(expr.dt$gene)
genes.to.plot <- c("Foxa2","Tfap2a","Mesp1")

for (i in genes.to.plot) {
  
  to.plot <- expr.dt[gene==i] %>% setorder(-expr.atlas)
  
  # Dynamic axes
  p <- ggplot(to.plot, aes(x=expr.atlas, y=expr.query)) +
    geom_point(aes(fill=celltype), color="black", shape=21, size=5, stroke=0.5) +
    stat_cor(method = "pearson") +
    # geom_abline(slope=1, intercept=0, linetype="dashed", color="black") +
    stat_smooth(method="lm", color="black", alpha=0.1, size=0.25) +
    scale_fill_manual(values=opts$celltype.colors) +
    # ggrepel::geom_text_repel(aes(label=celltype), segment.size=0.15, size=2.5, data=to.plot[expr.atlas>8 | expr.query>8]) +
    ggrepel::geom_text_repel(aes(label=gsub("_"," ",celltype)), segment.size=0.15, size=3.5, data=head(to.plot,n=7)) +
    labs(x="RNA expression (PijuanSala2019)", y="RNA expression (This study)", title=i) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color="black"),
      legend.position = "none"
    )

  pdf(file.path(io$outdir,sprintf("per_gene/%s_scatterplot_PijuanSala2019_comparison.pdf",i)), width = 6, height = 5)
  print(p)
  dev.off()
  
  # Fix axes
  # max.lim <- max(c(to.plot$expr.atlas,to.plot$expr.query))
  # min.lim <- 0 # min(c(to.plot$expr.atlas,to.plot$expr.query))
  # p <- ggplot(to.plot, aes(x=expr.atlas, y=expr.query)) +
  #   geom_point(aes(fill=celltype), color="black", shape=21, size=3, stroke=0.5) +
  #   geom_abline(slope=1, intercept=0) + 
  #   coord_cartesian(xlim=c(min.lim,max.lim), ylim=c(min.lim,max.lim)) +
  #   scale_fill_manual(values=opts$celltype.colors) +
  #   ggrepel::geom_text_repel(aes(label=celltype), segment.size=0.15, size=2.5, data=to.plot[expr.atlas>0.1 | expr.query>0.1]) +
  #   labs(x="RNA expression (Atlas)", y="RNA expression (This study)", title=i) +
  #   theme_classic() + 
  #   theme(
  #     plot.title = element_text(hjust = 0.5),
  #     axis.text = element_text(color="black"),
  #     legend.position = "none"
  #   )
  # 
  # pdf(sprintf("%s/%s_scatterplot_PijuanSala2019_comparison.pdf",io$outdir,i), width = 6, height = 5)
  # print(p)
  # dev.off()
}


################################
## Correlations per cell type ##
###############################

cor_celltype.dt <- expr.dt %>% 
  .[, cor.test(x=expr.atlas, y=expr.query)[c("estimate","p.value")], by = c("celltype")] %>%
  .[, padj_fdr := p.adjust(p.value, method="fdr")] %>%
  .[, sig := padj_fdr <= 0.1] %>%
  setorder(padj_fdr)

to.plot <- cor_celltype.dt %>% .[,celltype:=factor(celltype,levels=rev(opts$celltypes))]

p <- ggplot(to.plot, aes(x=celltype, y=estimate, fill=celltype)) +
  geom_bar(stat="identity", width=0.45, alpha=0.9) +
  geom_point(shape=21, size=3.5) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  coord_flip() +
  labs(x="", y="Correlation coefficient") +
  theme(
    axis.text = element_text(colour="black",size=rel(0.75)),
    axis.title.x = element_text(colour="black",size=rel(1.0)),
    legend.position="none"
  )

pdf(file.path(io$outdir,"PijuanSala2019_comparison_correlation_per_celltype.pdf"), width = 6, height = 7)
print(p)
dev.off()


################################
## Scatterplots per cell type ##
################################

# celltypes.to.plot <- unique(expr.dt$celltype)
celltypes.to.plot <- c("Notochord","Neural_crest","Nascent_mesoderm")

for (i in celltypes.to.plot) {
  
  to.plot <- expr.dt[celltype==i] %>% setorder(-expr.atlas)
  
  p <- ggplot(to.plot, aes(x=expr.atlas, y=expr.query)) +
    geom_point(color="black", fill="gray90", alpha=0.75, shape=21, size=1.5, stroke=0.5) +
    stat_cor(method = "pearson") +
    # geom_abline(slope=1, intercept=0, linetype="dashed", color="black") +
    stat_smooth(method="lm", color="black", alpha=0.15, size=0.5) +
    ggrepel::geom_text_repel(aes(label=gene), segment.size=0.15, size=3, data=head(to.plot,n=10)) +
    labs(x="RNA expression (PijuanSala2019)", y="RNA expression (This study)", title=i) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color="black", size=rel(0.8)),
      legend.position = "none"
    )
  
  pdf(file.path(io$outdir,sprintf("per_celltype/%s_scatterplot_PijuanSala2019_comparison.pdf",i)), width = 6, height = 5)
  print(p)
  dev.off()
}

