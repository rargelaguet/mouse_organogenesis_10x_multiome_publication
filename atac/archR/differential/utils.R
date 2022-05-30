gg_volcano_plot <- function(to.plot, top_genes=10, xlim=NULL, ylim=NULL, label_groups = NULL) {
  
  negative_hits <- to.plot[sig==TRUE & MeanDiff<0,idx]
  positive_hits <- to.plot[sig==TRUE & MeanDiff>0,idx]
  all <- nrow(to.plot)
  
  if (is.null(xlim))
    xlim <- max(abs(to.plot$MeanDiff), na.rm=T)
  if (is.null(ylim))
    ylim <- max(-log10(to.plot$FDR), na.rm=T)
  
  p <- ggplot(to.plot, aes(x=MeanDiff, y=-log10(FDR))) +
    # ggrastr::geom_point_rast(aes(color=sig), size=1) +
    geom_point(aes(color=sig), size=1) +
    # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
    geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange") +
    scale_color_manual(values=c("black","red")) +
    # scale_x_continuous(limits=c(-xlim-10,xlim+10)) +
    scale_x_continuous(limits=c(-xlim,xlim)) +
    scale_y_continuous(limits=c(0,ylim+2.5)) +
    labs(x="Accessibility difference (%)", y=expression(paste("-log"[10],"(FDR)"))) +
    annotate("text", x=0, y=ylim+1, size=7, label=sprintf("(%d)", all)) +
    annotate("text", x=-50, y=ylim+2, size=7, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=50, y=ylim+2, size=7, label=sprintf("%d (+)",length(positive_hits))) +
    # ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=top_genes), aes(x=MeanDiff, y=-log10(FDR), label=symbol), size=5) +
    theme_classic() +
    theme(
      axis.text = element_text(size=rel(1.00), color='black'),
      axis.title = element_text(size=rel(1.50), color='black'),
      # axis.title = element_text(),
      legend.position="none"
    )

  if (length(label_groups)>0) {
    p <- p +
      annotate("text", x=-70, y=0, size=4.5, label=sprintf("Up in %s",label_groups[2])) +
      annotate("text", x=70, y=0, size=4.5, label=sprintf("Up in %s",label_groups[1]))
  }

  return(p)
}
