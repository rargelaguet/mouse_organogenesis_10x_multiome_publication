source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$limma_results_sturm <- "/Users/rargelaguet/data/Ahangari2022_NHLF_tgfb/rna/results/limma/limma_output_passage.csv.gz"
io$limma_results_hdf <- "/Users/rargelaguet/data/hdf_reik_batch2/rna/results/limma/age/39_donors/limma_results_horvath2_clock.txt"
io$outdir <- file.path(io$basedir,"rna/results/limma/comparison_sturm_vs_hdf_reik"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options


########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

gene_sets <- list()
gene_sets[["Transcription Factors"]] <- fread(io$TFs, header=F)[[1]]
gene_sets[["Surface markers"]] <- unique(fread(io$cell_surface_markers, header = F)[[1]])
gene_sets[["Epigenetic factors"]] <- unique(fread(io$epigenetic_factors, header=F)[[1]])
gene_sets[["Ribosomal"]] <-  fread(io$ribo_genes, header=F)[[1]]
gene_sets[["Mitochondrial"]] <- fread(io$mito_genes, header=F)[[1]]
gene_sets[["DNA repair"]] <- fread(io$dna_repair_genes, header=F)[[1]]
gene_sets[["Cell cycle"]] <- fread(io$cell_cycle_g2m_genes, header=F)[[1]]
gene_sets[["ECM"]] <- fread(io$ecm_genes, header=F)[[1]]
# histones_genes <- fread(io$histones_genes, header=F)[[1]]


########################
## Load limma results ##
########################

limma_results_sturm.dt <- fread(io$limma_results_sturm) %>%
  .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
  .[,log_pval:=-log10(pval_adj+1e-100)] %>%
  .[,c("pval_adj","ens_id"):=NULL]

limma_results_hdf_reik.dt <- fread(io$limma_results_hdf) %>%
  .[is.na(coef),c("coef","pval_adj"):=list(0,1)] %>%
  .[,log_pval:=-log10(pval_adj+1e-100)] %>%
  .[,c("pval","pval_adj"):=NULL]


limma_results.dt <- merge(limma_results_sturm.dt, limma_results_hdf_reik.dt, by=c("feature"), suffixes = c("_A", "_B")) %>%
  .[,transcript_type:="Genes"] %>%
  .[feature%in%gene_sets[["Surface markers"]],transcript_type:="Surface markers"] %>%
  .[feature%in%gene_sets[["Epigenetic factors"]],transcript_type:="Epigenetic factors"] %>%
  .[feature%in%gene_sets[["Ribosomal"]],transcript_type:="Ribosomal"] %>%
  .[feature%in%gene_sets[["Mitochondrial"]],transcript_type:="Mitochondrial"] %>%
  .[feature%in%gene_sets[["DNA repair"]],transcript_type:="DNA repair"] %>%
  .[feature%in%gene_sets[["Cell cycle"]],transcript_type:="Cell cycle"] %>%
  .[feature%in%gene_sets[["ECM"]],transcript_type:="ECM"] %>%
  .[feature%in%gene_sets[["Transcription Factors"]],transcript_type:="Transcription Factors"] %>%
  .[,transcript_type:=factor(transcript_type, levels=c("Genes",names(gene_sets)))]

##########
## Plot ##
##########

to.plot <- limma_results.dt %>%
  .[,tmp:=abs(coef_A)*abs(coef_B)] %>% setorder(-tmp) %>% .[,tmp:=NULL] %>% 
  .[abs(coef_A)>0.01 | abs(coef_B)>0.003]

xlim.max <- max(abs(to.plot$coef_A))
ylim.max <- max(abs(to.plot$coef_B))

genes.to.highlight <- c("TXNIP")

p <- ggplot(to.plot, aes(x=coef_A, y=coef_B)) +
  geom_point(size=0.15, color="gray50", data=to.plot[!(abs(coef_B)>=0.01 | abs(coef_A)>=0.01)][!feature%in%genes.to.highlight]) +
  geom_point(size=0.45, color="gray50", data=to.plot[abs(coef_B)>=0.01 | abs(coef_A)>=0.01][!feature%in%genes.to.highlight]) +
  geom_point(size=0.45, color="purple", data=to.plot[feature%in%genes.to.highlight]) +
  # stat_smooth(method="lm") +
  stat_cor(aes(label = after_stat(r.label)),  method = "pearson") +
  geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
  geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
  labs(x="limma coefficient (Sturm2022)", y="limma coefficient (Reik)", title="") +
  geom_text_repel(data=to.plot[feature%in%genes.to.highlight], aes(x=coef_A, y=coef_B, label=feature), color="purple", size=2.50, segment.size=0.1, max.overlaps = Inf) +
  facet_wrap(~transcript_type, nrow=2) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5, size=rel(1.0)),
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

for (j in unique(to.plot$transcript_type)) {
  p <- p +
    geom_text_repel(data=head(to.plot[transcript_type==j & coef_B>=0.01 & coef_A>=0.01 & !feature%in%genes.to.highlight],n=15), aes(x=coef_A, y=coef_B, label=feature), color="black", size=2.50, segment.size=0.1, max.overlaps = Inf) +
    geom_text_repel(data=head(to.plot[transcript_type==j & coef_B<(-0.01) & coef_A<=(-0.01) & !feature%in%genes.to.highlight],n=15), aes(x=coef_A, y=coef_B, label=feature), color="black", size=2.50, segment.size=0.1, max.overlaps = Inf)
}

pdf(file.path(io$outdir,"limma_scatterplot_sturm_vs_reik.pdf"), width=14, height=10)
print(p)
dev.off()

#######################
## Plot per gene set ##
#######################

outdir <- file.path(io$outdir,"per_gene_set"); dir.create(outdir, showWarnings = F, recursive = T)

for (i in names(gene_sets)) {
  
  to.plot.i <- to.plot[transcript_type==i]
  
  p <- ggplot(to.plot.i, aes(x=coef_A, y=coef_B)) +
    geom_point(size=0.15, color="gray50", data=to.plot.i[!(abs(coef_B)>=0.01 | abs(coef_A)>=0.01)][!feature%in%genes.to.highlight]) +
    geom_point(size=0.45, color="gray50", data=to.plot.i[abs(coef_B)>=0.01 | abs(coef_A)>=0.01][!feature%in%genes.to.highlight]) +
    geom_point(size=0.45, color="purple", data=to.plot.i[feature%in%genes.to.highlight]) +
    # stat_smooth(method="lm") +
    stat_cor(aes(label = after_stat(r.label)),  method = "pearson") +
    geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
    geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
    geom_text_repel(data=head(to.plot.i[transcript_type==i & coef_B>=0.01 & coef_A>=0.01 & !feature%in%genes.to.highlight],n=15), aes(x=coef_A, y=coef_B, label=feature), color="black", size=2.50, segment.size=0.1, max.overlaps = Inf) +
    geom_text_repel(data=head(to.plot.i[transcript_type==i & coef_B<(-0.01) & coef_A<=(-0.01) & !feature%in%genes.to.highlight],n=15), aes(x=coef_A, y=coef_B, label=feature), color="black", size=2.50, segment.size=0.1, max.overlaps = Inf) +
    geom_text_repel(data=head(to.plot.i[transcript_type==i & abs(coef_B)<0.01 & coef_A>=(0.01) & !feature%in%genes.to.highlight] %>% setorder(-coef_A),n=15), aes(x=coef_A, y=coef_B, label=feature), color="black", size=2.50, segment.size=0.1, max.overlaps = Inf) +
    geom_text_repel(data=head(to.plot.i[transcript_type==i & coef_B>(0.01) & abs(coef_A)<0.01 & !feature%in%genes.to.highlight] %>% setorder(-coef_B),n=15), aes(x=coef_A, y=coef_B, label=feature), color="black", size=2.50, segment.size=0.1, max.overlaps = Inf) +
    labs(x="limma coefficient (Sturm2022)", y="limma coefficient (Reik)", title="") +
    geom_text_repel(data=to.plot.i[feature%in%genes.to.highlight], aes(x=coef_A, y=coef_B, label=feature), color="purple", size=2.50, segment.size=0.1, max.overlaps = Inf) +
    facet_wrap(~transcript_type, nrow=2) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust=0.5, size=rel(1.0)),
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black'),
      legend.position="none"
    )
  
  pdf(file.path(outdir,sprintf("limma_scatterplot_sturm_vs_reik_%s.pdf",gsub(" ","-",i))), width=8, height=7)
  print(p)
  dev.off()
}
