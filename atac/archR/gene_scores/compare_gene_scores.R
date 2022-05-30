#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

io$pseudobulk.dir <- "/Users/ricard/data/gastrulation_multiome_10x/processed/atac/archR/pseudobulk"

###############
## Load data ##
###############

gene_scores_v1.se <- readRDS(paste0(io$pseudobulk.dir,"/pseudobulk_GeneScoreMatrix_nodistal_summarized_experiment.rds"))
gene_scores_v2.se <- readRDS(paste0(io$pseudobulk.dir,"/pseudobulk_GeneScoreMatrix_nodistal_v2_normalised_summarized_experiment.rds"))

rownames(gene_scores_v1.se) <- rowData(gene_scores_v1.se)$name
rownames(gene_scores_v2.se) <- rowData(gene_scores_v2.se)$name

# Convert to data.table 
gene_scores_v1.dt <- assay(gene_scores_v1.se) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="value") %>%
  .[,class:="v1"]


gene_scores_v2.dt <- assay(gene_scores_v2.se) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="value") %>%
  .[,class:="v2"]

gene_scores.dt <- rbind(gene_scores_v1.dt,gene_scores_v2.dt) %>%
	dcast(gene+celltype~class)

#################
## Correlation ##
#################

cor.dt <- gene_scores.dt %>% copy %>%
  .[, .(V1 = unlist(cor.test(v1,v2)[c("estimate", "p.value")])), by = c("gene")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene ~ para, value.var = "V1") %>%
  setorder(p, na.last = T)

# Save
# fwrite(cor.dt, paste0(io$outdir,"/per_gene/cor_rna_vs_acc_pseudobulk.txt.gz"), sep="\t", quote=F)


##################
## Volcano plot ##
##################

to.plot <- cor.dt

negative_hits <- to.plot[sig==TRUE & r<0,gene]
positive_hits <- to.plot[sig==TRUE & r>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$r), na.rm=T)
ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)

p <- ggplot(to.plot, aes(x=r, y=-log10(padj_fdr+1e-100))) +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.25) +
  ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & r>0],n=50),
                           aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3, max.overlaps=100) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & r<0],n=10),
                           aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3, max.overlaps=100) +
  scale_color_manual(values=c("black","red")) +
  scale_size_manual(values=c(0.75,1.25)) +
  scale_x_continuous(limits=c(-xlim-0.2,xlim+0.2)) +
  scale_y_continuous(limits=c(0,ylim+6)) +
  annotate("text", x=0, y=ylim+6, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.15, y=ylim+6, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.15, y=ylim+6, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

##########
## Plot ##
##########

to.plot <- assay(GeneScoreMatrix_pseudobulk.se[selected.genes,]) %>% as.matrix %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype") %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

ggbarplot(to.plot, x="celltype", y="value", fill="celltype") +
  facet_wrap(~gene) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Accessibility") +
  geom_hline(yintercept=0, colour="black", linetype="dashed") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    # axis.text.x = element_text(size=rel(0.75)),
    axis.text.x = element_blank(),
    # axis.title.y = element_text(size=rel(0.75)),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  )
