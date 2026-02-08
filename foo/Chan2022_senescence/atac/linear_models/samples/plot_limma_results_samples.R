# here::i_am("differential/differential_gene_expr.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# io$limma_results <- file.path(io$basedir,"results/atac/limma/all_samples/limma_results.rds")
# io$outdir <- file.path(io$basedir,"results/atac/limma/all_samples/pdf"); dir.create(io$outdir, showWarnings = F)
io$limma_results <- file.path(io$basedir,"results/atac/limma/control_samples/limma_results.rds")
io$outdir <- file.path(io$basedir,"results/atac/limma/control_samples/pdf"); dir.create(io$outdir, showWarnings = F)

########################
## Load limma results ##
########################

limma.list <- readRDS(io$limma_results)

###################
## Volcano plots ##
###################

constrasts <- names(limma.list)

# i <- constrasts[1]
for (i in constrasts) {
  
  # groupA <- "mix"; groupB <- "nr1"
  
  to.plot <- limma.list[[i]] %>%
    .[,sig:=pval_adj<0.01 & abs(logFC)>=0.1] %>%
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
    scale_size_manual(values=c(0.20,0.60)) +
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
  
  pdf(file.path(io$outdir,sprintf("%s_atac_volcano_plot_Chan2022_senescence.pdf",i)), width=10, height=7)
  print(p)
  dev.off()
}

##################################
## Plot directionality of logFC ##
##################################

# to.plot <- diff_hgps.dt %>% 
#   .[family%in%opts$families.to.plot] %>%
#   .[,sign:=ifelse(logFC>0,"Downregulated in L1-ASO treated cells","Upregulated in L1-ASO treated cells")] %>%
#   .[,.(N=sum(sig==T)),by=c("class","family","sign")]
# 
# p <- ggplot(to.plot, aes(x=family, y=N, fill=sign)) +
#   geom_bar(stat="identity", color="black") +
#   scale_fill_brewer(palette = "Dark2") +
#   facet_grid(~class, scales = "free_x", space = "free_x") +
#   labs(x="", y="Number of differentially expressed TEs") +
#   guides(x = guide_axis(angle = 90)) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(size=rel(1.0), color="black"),
#     axis.text.y = element_text(size=rel(1.0), color="black"),
#     legend.position = "top",
#     legend.title = element_blank()
#   )
# 
# pdf(file.path(io$outdir,"diff_directionality.pdf"), width=10, height=6)
# print(p)
# dev.off()

##############################
## Mean expression vs logFC ##
##############################

# # to.plot <- diff_hgps.dt %>% copy %>% .[,mean:=(mean_groupA+mean_groupB)/2] %>% .[mean<=15]
# to.plot <- diff_wrn.dt %>% copy %>% .[,mean:=(mean_groupA+mean_groupB)/2] %>% .[mean<=15]
# 
# to.plot[,highlight:=abs(logFC)>=0.1]
# xlim <- max(abs(to.plot$mean))
# ylim <- max(2,max(abs(to.plot$logFC)))
# 
# p <- ggplot(to.plot, aes(x=mean, y=logFC)) +
#   geom_jitter(aes(size=highlight), width=0.01, height=0.01) +
#   scale_size_manual(values=c(0.25,1.50)) +
#   geom_segment(x=0, xend=xlim, y=0, yend=0, color="orange", size=0.1) +
#   ggrepel::geom_text_repel(data=to.plot[abs(logFC)>0.25], aes(x=mean, y=logFC, label=feature), size=2.2, segment.size=0.1, max.overlaps = Inf) +
#   labs(x="Mean expression (WRN)", y="logFC (WRN)") +
#   annotate("text", x=1, y=ylim/1.1, size=4, label="Down in L1-ASO") +
#   annotate("text", x=1, y=-ylim/1.1, size=4, label="Up in L1-ASO") +
#   coord_cartesian(ylim=c(-ylim,ylim)) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(size=rel(1.0), color="black"),
#     axis.text.y = element_text(size=rel(1.0), color="black"),
#     legend.position = "none",
#     legend.title = element_blank()
#   )
# 
# pdf(file.path(io$outdir,"WRN_mean_vs_logFC_TE_expr_DellaValle2022.pdf"), width=9, height=7)
# print(p)
# dev.off()
