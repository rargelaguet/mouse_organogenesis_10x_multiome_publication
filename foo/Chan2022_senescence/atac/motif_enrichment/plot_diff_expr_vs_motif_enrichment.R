source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

# I/O
io$feature_metadata <- file.path(io$basedir,"processed/atac/archR/peak_metadata.tsv.gz")
io$limma_motif_enrichment <- file.path(io$basedir,sprintf("results/atac/archR/motif_enrichment/%s/motif_enrichment_limma_acc.txt.gz",opts$motif_annotation))
io$motif2gene <- file.path(io$basedir,sprintf("processed/atac/bigwig/deeptools/ATAC_peaks/motif_annotation/%s/%s_motif2gene.txt.gz",opts$motif_annotation,opts$motif_annotation))
io$limma_expr <- file.path(io$basedir,"results/bulk_rna/limma/limma_results_horvath2_clock.txt")
io$outdir <- file.path(io$basedir,"results/rna_atac/motif_enrichment_vs_diff_expr")

###########################
## Load feature metadata ##
###########################

feature_metadata.dt <- fread(io$feature_metadata)

##########################
## Load expression data ##
##########################

# sce <- readRDS(io$sce)

############################
## Load motif annotations ##
############################

motif2gene.dt <- fread(io$motif2gene) %>%
  .[,c("motif","gene")]


##########################################
## Load differential expression results ##
##########################################

limma_expr.dt <- fread(io$limma_expr) %>% .[feature%in%TFs]

###################################
## Load motif enrichment results ##
###################################

motif_enrichment.dt <- fread(io$limma_motif_enrichment) %>% 
  # .[contrast=="horvath2_clock"] %>%
  .[contrast=="chronological_age"] %>%
  merge(motif2gene.dt, by="motif") %>%
  .[,log_pval:=-log10(pval)] %>% .[log_pval>=150,log_pval:=150]

######################################################
## Plot differential expression vs motif enrichment ##
######################################################

to.plot <- merge(
  limma_expr.dt[,c("feature","coef")] %>% setnames(c("gene","logFC_rna")),
  motif_enrichment.dt[sign=="Open",c("gene","log_pval")] %>% setnames(c("gene","log_pval_acc")),
  by = c("gene")
) %>% .[,tmp:=logFC_rna*log_pval_acc]

p1 <- ggplot(to.plot, aes(x=log_pval_acc, y=logFC_rna)) +
  # geom_point(aes(size=tmp)) +
  geom_point(size=1) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  # scale_size_continuous(range = c(0.5,1.5)) +
  labs(y="Log fold change (expression)", x=expression(paste("-log"[10],"(p.value) motif enrichment")), title="ATAC peaks open in old donors") +
  # annotate("text", x=20, y=2, size=3, label=sprintf("Up in %s",j)) +
  ggrepel::geom_text_repel(data=head(to.plot[abs(logFC_rna)>=0.030],n=50), aes(x=log_pval_acc, y=logFC_rna, label=gene), segment.size=0.05, size=2.25, max.overlaps = Inf) +
  ggrepel::geom_text_repel(data=head(to.plot[log_pval_acc>=20],n=50), aes(x=log_pval_acc, y=logFC_rna, label=gene), segment.size=0.05, size=2.25, max.overlaps = Inf) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5, size=rel(1.0)),
    plot.subtitle = element_text(hjust=0.5, size=rel(0.75)),
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position = "none"
  )

to.plot <- merge(
  limma_expr.dt[,c("feature","coef")] %>% setnames(c("gene","logFC_rna")),
  motif_enrichment.dt[sign=="Closed",c("gene","log_pval")] %>% setnames(c("gene","log_pval_acc")),
  by = c("gene")
) %>% .[,tmp:=logFC_rna*log_pval_acc]

p2 <- ggplot(to.plot, aes(x=log_pval_acc, y=logFC_rna)) +
  geom_point(size=1) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  # scale_size_continuous(range = c(0.5,1.5)) +
  labs(y="Log fold change (expression)", x=expression(paste("-log"[10],"(p.value) motif enrichment")), title="ATAC peaks closed in old donors") +
  # annotate("text", x=20, y=2, size=3, label=sprintf("Up in %s",j)) +
  ggrepel::geom_text_repel(data=head(to.plot[abs(logFC_rna)>=0.030],n=50), aes(x=log_pval_acc, y=logFC_rna, label=gene), segment.size=0.05, size=2.25, max.overlaps = Inf) +
  ggrepel::geom_text_repel(data=head(to.plot[log_pval_acc>=20],n=50), aes(x=log_pval_acc, y=logFC_rna, label=gene), segment.size=0.05, size=2.25, max.overlaps = Inf) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5, size=rel(1.0)),
    plot.subtitle = element_text(hjust=0.5, size=rel(0.75)),
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position = "none"
  )

p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1) 

pdf(file.path(io$outdir,"limma_expr_vs_motif_enrichment_chronological_age_fibros.pdf"), width=12, height=6)
print(p)
dev.off()
