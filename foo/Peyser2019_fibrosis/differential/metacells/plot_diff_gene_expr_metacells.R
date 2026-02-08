source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$group_variable <- "condition"
opts$celltypes <- c("Type_II_alveolar", "dendritic_cells", "monocytes", "Type_I_alveolar", "mesothelium", "endothelial", "macrophages", "B_cells", "fibroblasts")
# opts$celltypes <- c("fibroblasts")

# I/O
io$diff_results <- file.path(io$basedir,sprintf("results/differential/metacells/%s/parsed/diff_expr_results.txt.gz",opts$group_variable))
io$outdir <- file.path(io$basedir,"results/differential/metacells/group_manual/pdf"); dir.create(io$outdir, showWarnings = F)
io$outdir.volcano <- file.path(io$basedir,sprintf("results/differential/metacells/%s/pdf/volcano_plots",opts$group_variable)); dir.create(io$outdir.volcano, showWarnings = F, recursive = T)
io$outdir.scatter <- file.path(io$basedir,sprintf("results/differential/metacells/%s/pdf/scatter_plots",opts$group_variable)); dir.create(io$outdir.scatter, showWarnings = F, recursive = T)
# io$outdir.xing2020 <- file.path(io$basedir,sprintf("results/differential/metacells/%s/pdf/xing2020_vs_bj",opts$group_variable)); dir.create(io$outdir.xing2020, showWarnings = F, recursive = T)

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

gene_sets <- list()
gene_sets[["Transcription Factors"]] <- fread(io$TFs, header=F)[[1]] %>% str_to_title(.)
gene_sets[["Surface markers"]] <- unique(fread(io$cell_surface_markers, header = F)[[1]]) %>% str_to_title(.)
gene_sets[["Epigenetic factors"]] <- unique(fread(io$epigenetic_factors, header=F)[[1]]) %>% str_to_title(.)
gene_sets[["Ribosomal"]] <-  fread(io$ribo_genes, header=F)[[1]] %>% str_to_title(.)
gene_sets[["Mitochondrial"]] <- fread(io$mito_genes, header=F)[[1]]
# gene_sets[["DNA repair"]] <- fread(io$dna_repair_genes, header=F)[[1]] %>% str_to_title(.)
gene_sets[["Cell cycle"]] <- fread(io$cell_cycle_g2m_genes, header=F)[[1]] %>% str_to_title(.)
gene_sets[["ECM"]] <- fread(io$ecm_genes, header=F)[[1]] %>% str_to_title(.)
# histones_genes <- fread(io$histones_genes, header=F)[[1]] %>% str_to_title(.)

##########################################
## Load differential expression results ##
##########################################

diff.dt <- fread(io$diff_results) %>%
  setnames("gene","feature") %>%
  # .[groupA%in%opts$groups & groupB%in%opts$groups] %>%
  .[is.na(logFC),c("logFC","padj_fdr"):=list(0,1)] %>%
  .[,log_pval:=-log10(padj_fdr+1e-300)] %>%
  .[,transcript_type:="Genes"] %>%
  .[feature%in%gene_sets[["Surface markers"]],transcript_type:="Surface markers"] %>%
  .[feature%in%gene_sets[["Epigenetic factors"]],transcript_type:="Epigenetic factors"] %>%
  .[feature%in%gene_sets[["Ribosomal"]],transcript_type:="Ribosomal"] %>%
  .[feature%in%gene_sets[["Mitochondrial"]],transcript_type:="Mitochondrial"] %>%
  # .[feature%in%gene_sets[["DNA repair"]],transcript_type:="DNA repair"] %>%
  .[feature%in%gene_sets[["Cell cycle"]],transcript_type:="Cell cycle"] %>%
  .[feature%in%gene_sets[["ECM"]],transcript_type:="ECM"] %>%
  .[feature%in%gene_sets[["Transcription Factors"]],transcript_type:="Transcription Factors"] %>%
  .[,transcript_type:=factor(transcript_type, levels=c("Genes",names(gene_sets)))]

##########################
## Filter diff results ##
##########################

diff.dt <- diff.dt[!feature%in%genes.to.remove]

# diff.dt[feature=="Txnip" & celltype=="fibroblasts"]

###################
## Volcano plots ##
###################

i <- "Saline_Control"; j <- "Bleomycin_day11"
# i <- "Saline_Control"; j <- "Bleomycin_nintedanib_day11"
# i <- "Bleomycin_day11"; j <- "Bleomycin_nintedanib_day11"
for (c in opts$celltypes) {
  
  diff_filt.dt <- diff.dt[groupA==i & groupB==j & celltype==c] %>% 
    .[,sig:=padj_fdr<0.01 & abs(logFC)>=0.50]
  
  tmp <- diff_filt.dt[,.(N=as.character(.N), Npos=as.character(sum(sig==TRUE & logFC>0)), Nneg=as.character(sum(sig==TRUE & logFC<0))),by="transcript_type"]
  tmp %>% .[,Npos:=sprintf("%s (+)",Npos)] %>% .[,Nneg:=sprintf("%s (-)",Nneg)]
  
  to.plot <- diff_filt.dt %>% 
    .[!grepl("Rik",feature)] %>%
    .[abs(logFC)>0.05 & padj_fdr<1] %>%
    .[,log_pval:=-log10(padj_fdr+1e-300)]
  
  negative_hits <- to.plot[sig==TRUE & logFC<0,feature]
  positive_hits <- to.plot[sig==TRUE & logFC>0,feature]
  all <- nrow(to.plot)
  
  xlim <- max(c(2,abs(to.plot$logFC)), na.rm=T)
  ylim <- max(c(10,-log10(to.plot$padj_fdr+1e-100)), na.rm=T)
  
  p <- ggplot(to.plot, aes(x=logFC, y=log_pval)) +
    geom_point(aes(size=sig), color="gray40") +
    geom_segment(x=0, xend=0, y=0, yend=ylim, color="orange", linewidth=0.15) +
    # scale_color_manual(values=c("black","red")) +
    scale_size_manual(values=c(0.25,0.50)) +
    scale_x_continuous(limits=c(-xlim-0.5,xlim+0.5)) +
    scale_y_continuous(limits=c(0,ylim+2)) +
    labs(x="Log fold change", y=expression(paste("-log"[10],"(p.value)")), title=sprintf("%s vs %s: %s",i,j,c)) +
    geom_text(aes(label=N), x=0, y=ylim+0.5, size=3, data=tmp) +
    geom_text(aes(label=Nneg), x=-xlim+(xlim/3), y=ylim+0.5, size=3, data=tmp) +
    geom_text(aes(label=Npos), x=xlim-(xlim/3), y=ylim+0.5, size=3, data=tmp) +
    # annotate("text", x=-xlim+(xlim/5), y=0, size=2, label=sprintf("Up in %s",i)) +
    # annotate("text", x=xlim-(xlim/5), y=0, size=2, label=sprintf("Up in %s",j)) +
    facet_wrap(~transcript_type, nrow=2) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust=0.5, size=rel(1.0)),
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black'),
      legend.position="none"
    )
  
  for (k in unique(to.plot$transcript_type)) {
    p <- p + 
      geom_text_repel(data=head(to.plot[transcript_type==k & sig==TRUE & logFC>0],n=10), aes(x=logFC, y=log_pval, label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[transcript_type==k & sig==TRUE & logFC<0],n=10), aes(x=logFC, y=log_pval, label=feature), segment.size=0.05, size=2.50, max.overlaps = Inf)
  }
  
  pdf(file.path(io$outdir.volcano,sprintf("volcano_plot_%s_%s_vs_%s_Peyser2019.pdf",c,i,j)), width=10, height=8)
  print(p)
  dev.off()
}


##################
## Scatterplots ##
##################

genes.to.highlight <- c()

for (c in opts$celltypes) {
  
  diff_1.dt <- diff.dt[groupA=="Saline_Control" & groupB=="Bleomycin_day11" & celltype==c,c("feature","transcript_type","logFC")]
  diff_2.dt <- diff.dt[groupA=="Saline_Control" & groupB=="Bleomycin_nintedanib_day11" & celltype==c,c("feature","transcript_type","logFC")]
  
  to.plot <- merge(diff_1.dt, diff_2.dt, by=c("feature","transcript_type"), suffixes=c("_A","_B")) %>%
    .[abs(logFC_A)>0.01 & abs(logFC_B)>0.01] %>%
    .[,tmp:=abs(logFC_A)] %>% setorder(-tmp) %>% .[,tmp:=NULL]
  
  # to.plot <- to.plot[transcript_type%in%c("Surface markers")]
  
  ylim.max <- max(abs(to.plot$logFC_A))
  xlim.max <- max(abs(to.plot$logFC_B))
  
  p <- ggplot(to.plot, aes(x=logFC_A, y=logFC_B)) +
    stat_cor(aes(label = after_stat(r.label)),  method = "pearson") +
    geom_point(size=0.15, color="gray30") +
    # geom_jitter(size=0.45, width=0.01, height=0.01, color="gray30", data=to.plot[abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5]) +
    # geom_jitter(size=0.15, width=0.01, height=0.01, color="gray30", data=to.plot[!(abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5)][!feature%in%genes.to.highlight]) +
    # geom_jitter(size=0.45, width=0.01, height=0.01, color="gray30", data=to.plot[abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5][!feature%in%genes.to.highlight]) +
    geom_point(size=0.45, color="purple", data=to.plot[feature%in%genes.to.highlight]) +
    geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
    geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
    ggrepel::geom_text_repel(data=to.plot[feature%in%genes.to.highlight], aes(x=logFC_A, y=logFC_B, label=feature), color="purple", size=2.50, segment.size=0.1, max.overlaps = Inf) +
    labs(x="logFC (Saline vs Bleomycin_day11)", y="logFC (Saline vs Bleomycin_nintedanib_day11)", title="") +
    facet_wrap(~transcript_type, nrow=3) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size=rel(1.0), color="black"),
      axis.text.y = element_text(size=rel(1.0), color="black"),
      legend.position = "none",
      legend.title = element_blank()
    )
  
  for (k in unique(to.plot$transcript_type)) {
    p <- p +
      geom_text_repel(data=head(to.plot[!feature%in%genes.to.highlight][transcript_type==k & logFC_B>=0.01 & logFC_A>=1],n=7), aes(x=logFC_A, y=logFC_B, label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[!feature%in%genes.to.highlight][transcript_type==k & abs(logFC_A)<(1) & logFC_B>=0.01] %>% setorder(-logFC_B),n=7), aes(x=logFC_A, y=logFC_B, label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[!feature%in%genes.to.highlight][transcript_type==k & logFC_A>1 & abs(logFC_B)<=0.01] %>% setorder(-logFC_A),n=7), aes(x=logFC_A, y=logFC_B, label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[!feature%in%genes.to.highlight][transcript_type==k & logFC_A>1 & logFC_B<=(-0.01)],n=7), aes(x=logFC_A, y=logFC_B, label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[!feature%in%genes.to.highlight][transcript_type==k & abs(logFC_A)<(1) & logFC_B<=(-0.01)] %>% setorder(logFC_B),n=7), aes(x=logFC_A, y=logFC_B, label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[!feature%in%genes.to.highlight][transcript_type==k & logFC_B<(-0.01) & logFC_A<=(-1)],n=7), aes(x=logFC_A, y=logFC_B, label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[!feature%in%genes.to.highlight][transcript_type==k & logFC_A<(-1) & abs(logFC_B)<=(0.01)] %>% setorder(logFC_A),n=7), aes(x=logFC_A, y=logFC_B, label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf) +
      geom_text_repel(data=head(to.plot[!feature%in%genes.to.highlight][transcript_type==k & logFC_A<(-1) & logFC_B>0.01],n=7) %>% setorder(-logFC_A), aes(x=logFC_A, y=logFC_B, label=feature), segment.size=0.05, size=2.25, max.overlaps = Inf)
  }
  
  pdf(file.path(io$outdir.scatter,sprintf("scatterplot_differential_control_vs_bleomycine_%s_Peyser2019.pdf",c)), width=10, height=9)
  print(p)
  dev.off()
}
