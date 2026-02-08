# here::i_am("differential/differential_gene_expr.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$limma_results <- file.path(io$basedir,"results/proteomics/limma/limma_results.rds")
io$outdir <- file.path(io$basedir,"results/proteomics/limma/pdf"); dir.create(io$outdir, showWarnings = F)

########################
## Load limma results ##
########################

limma.list <- readRDS(io$limma_results)

###################
## Volcano plots ##
###################

contrasts <- names(limma.list)

# i <- contrasts[1]
for (i in contrasts) {
  
  to.plot <- limma.list[[i]] %>%
    .[,sig:=pval_adj<0.01 & abs(logFC)>=0.5] %>%
    setorder(pval_adj)
  
  negative_hits <- to.plot[sig==TRUE & logFC<0,feature]
  positive_hits <- to.plot[sig==TRUE & logFC>0,feature]
  all <- nrow(to.plot)
  
  xlim <- max(c(0.5,abs(to.plot$logFC)), na.rm=T)
  ylim <- max(c(10,-log10(to.plot$pval_adj+1e-300)), na.rm=T)
  
  p <- ggplot(to.plot, aes(x=logFC, y=-log10(pval_adj+1e-300))) +
    geom_point(aes(color=sig, size=sig)) +
    geom_segment(x=0, xend=0, y=0, yend=ylim, color="orange", size=0.5) +
    scale_color_manual(values=c("black","red")) +
    scale_size_manual(values=c(0.20,1)) +
    # scale_x_continuous(limits=c(-xlim,xlim+1.5)) +
    # scale_y_continuous(limits=c(0,ylim+2)) +
    labs(x="Log fold change", y=expression(paste("-log"[10],"(p.value)")), title=i) +
    annotate("text", x=0, y=ylim+1, size=4, label=sprintf("(%d)", all)) +
    annotate("text", x=-xlim, y=ylim+1, size=4, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=xlim, y=ylim+1, size=4, label=sprintf("%d (+)",length(positive_hits))) +
    # annotate("text", x=-xlim, y=0, size=3, label=sprintf("Up in %s",groupB)) +
    # annotate("text", x=xlim, y=0, size=3, label=sprintf("Up in %s",groupA)) +
    ggrepel::geom_text_repel(data=head(to.plot,n=50), aes(x=logFC, y=-log10(pval_adj+1e-300), label=feature), size=2.5, max.overlaps = Inf) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust=0.5),
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black'),
      legend.position="none"
    )
  
  pdf(file.path(io$outdir,sprintf("%s_proteomics_volcano_plot_Chan2022_senescence.pdf",i)), width=10, height=7)
  print(p)
  dev.off()
}
