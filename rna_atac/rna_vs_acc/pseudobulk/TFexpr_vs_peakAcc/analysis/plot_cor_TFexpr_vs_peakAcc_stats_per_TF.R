
#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/per_TF"); dir.create(io$outdir, showWarnings = F)

# Options
opts$motif_annotation <- "Motif_cisbp"

###############
## Load data ##
###############

io$tf2peak_cor.se <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
# io$tf2peak_cor.se <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment_lenient.rds")
tf2peak_cor.se <- readRDS(io$tf2peak_cor.se)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
# 
# if (grepl("ricard",Sys.info()['nodename'])) {
#   source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
# } else if (grepl("ebi",Sys.info()['nodename'])) {
#   source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
# } else {
#   stop("Computer not recognised")
# }

######################
## Load motifmatchR ##
######################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
} else {
  stop("Computer not recognised")
}


##############################################
## Load RNA vs chromRNA correlations per TF ##
##############################################

cor.dt <- fread(paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromRNA_pseudobulk.txt.gz")) %>%
# cor.dt <- fread(paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromVAR_pseudobulk.txt.gz")) %>%
  .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])]


###############################################################
## Plot histogram of the nº of peaks per TF before and after ##
###############################################################

# (TO-DO: USE MINIMUM SCORE THRESHOLD)
opts$min.corr <- 0.50

before_cor_filtering.dt <- data.table(
  TF = colnames(motifmatcher.se),
  N = colSums(assay(motifmatcher.se))
)

after_cor_filtering.dt <- data.table(
  TF = colnames(tf2peak_cor.se),
  # N = rowSums(dropNA2matrix(assay(tf2peak_cor.se,"pvalue"))<0.01,na.rm=T)
  # N = rowSums(abs(dropNA2matrix(assay(tf2peak_cor.se,"cor")))>=opts$min.corr,na.rm=T)
  N = Matrix::colSums(abs(assay(tf2peak_cor.se,"cor"))>=opts$min.corr,na.rm=T)
) 

to.plot <-  merge(before_cor_filtering.dt, after_cor_filtering.dt, by="TF", suffixes=c("_chromVAR","_chromRNA")) %>%
  melt(id.vars="TF") %>%
  .[,log_value:=log10(value+1)]

p <- ggdensity(to.plot, x="log_value", fill="variable", y="..density..") +
  geom_vline(aes(xintercept=V1), color="black", data=to.plot[,median(log_value),by="variable"], linetype="dashed") +
  labs(x="Number of peaks per TF motif (log10)", y="Density") +
  scale_fill_brewer(palette="Dark2") +
  scale_x_continuous(limits = c(0.5,5)) +
  # guides(color=F) + 
  theme(
    axis.title = element_text(size=rel(1)),
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )


pdf(sprintf("%s/density_number_peaks_per_TF_before_and_after.pdf",io$outdir), width=5, height=4)
print(p)
dev.off()

#################
## Arrow plots ##
#################

tmp <- data.table(
  TF = colnames(motifmatcher.se),
  N = colSums(assay(motifmatcher.se))
)

to.plot <- data.table(
  TF = colnames(tf2peak_cor.se),
  number_cor_peaks = colSums(dropNA2matrix(assay(tf2peak_cor.se,"pvalue"))<0.01,na.rm=T)
) %>% merge(tmp,by="TF")
  
TF.order <- tmp %>% setorder(-N) %>% head(n=50) %>% .$TF
to.plot2 <- to.plot %>% .[TF%in%TF.order] %>% melt(id.vars="TF") %>% .[,TF:=factor(TF,levels=TF.order)]

p <- ggscatter(to.plot2, x="TF", y="value", color="variable", size=1.5) +
  geom_line(arrow = arrow(length=unit(0.12,"cm"), ends="first", type = "closed")) +
  labs(x="", y="Number of correlation events with peaks") +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size=rel(0.75), color="black"),
    axis.text.x = element_text(size=rel(0.65), color="black"),
    axis.ticks.x = element_blank()
  )


pdf(sprintf("%s/arrow_number_correlation_events_per_TF.pdf",io$outdir), width = 8, height = 4)
print(p)
dev.off()

##########################################################################
## Plot number of correlation events per TF, split by peakType and sign ##
##########################################################################

tf2peak_cor.dt <- data.table(
  TF = colnames(tf2peak_cor.se),
  # number_no_cor = Matrix::colSums(assay(tf2peak_cor.se,"pvalue")>0.1,na.rm=T),
  number_positive_cor = Matrix::colSums(assay(tf2peak_cor.se,"cor")>0.30 & assay(tf2peak_cor.se,"pvalue")<0.10,na.rm=T),
  number_negative_cor = Matrix::colSums(assay(tf2peak_cor.se,"cor")<(-0.30) & assay(tf2peak_cor.se,"pvalue")<0.10,na.rm=T)
) 

tf2peak_cor.dt[,fraction_positive_cor:=number_positive_cor/(number_positive_cor+number_negative_cor)]

to.plot <- tf2peak_cor.dt[,c("TF","number_positive_cor","number_negative_cor")] %>% 
  melt(id.vars="TF", variable.name="cor_sign", value.name="N") %>%
  merge(cor.dt[,c("gene","cor_sign")] %>% setnames(c("TF","TF_type")))

TF.activators <- to.plot[TF_type=="Activator" & cor_sign=="number_positive_cor"] %>% setorder(N) %>% tail(n=25) %>% .$TF
TF.repressors <- to.plot[TF_type=="Repressor" & cor_sign=="number_negative_cor"] %>% setorder(N) %>% tail(n=25) %>% .$TF
# TF.mixed <- tf2peak_cor.dt[fraction_positive_cor>0.30 & fraction_positive_cor<0.70 & number_positive_cor>25] %>% setorder(number_positive_cor) %>% tail(n=25) %>% .$TF
TFs.to.plot <- c(TF.activators, TF.repressors)
to.plot.subset <- to.plot[TF%in%TFs.to.plot] %>% .[,TF:=factor(TF,levels=rev(TFs.to.plot))]
# to.plot.subset[TF%in%TF.mixed,TF_type:="Mixed"]

p <- ggbarplot(to.plot.subset, x="TF", y="N", fill="cor_sign") +
  labs(x="", y="Number of correlation events with ATAC peaks") +
  coord_flip() +
  facet_wrap(~TF_type, scales = "free") +
  theme_classic() +
  theme(
    strip.background =element_rect(fill=alpha("#DBDBDB", 0.65)),
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )


pdf(sprintf("%s/barplot_number_correlation_events_perTFType.pdf",io$outdir), width = 6, height = 5)
print(p)
dev.off()

##############################################
## Plot number of correlation events per TF ##
##############################################

# to.plot1 <- cor_dt[,.N,by="TF"] %>% 
#   setorder(-N) %>% 
#   .[,TF:=factor(TF,levels=TF)]
# 
# p <- ggplot(to.plot1 %>% head(n=100), aes_string(x="TF", y="N"), fill="gray70") +
#   geom_point(size=2) +
#   geom_segment(aes_string(xend="TF"), size=0.5, yend=0) +
#   coord_flip() +
#   labs(x="", y="Number of correlation events with peaks") +
#   theme_classic() +
#   theme(
#     axis.ticks.y = element_blank(),
#     axis.text = element_text(size=rel(0.75), color="black")
#   )
# 
# 
# pdf(sprintf("%s/barplot_number_correlation_events.pdf",io$outdir), width = 6, height = 10)
# print(p)
# dev.off()

##############################################
## Distribution of correlation coefficients ##
##############################################

pvalue.mtx <- dropNA2matrix(assay(tf2peak_cor.se,"pvalue"))
cor.mtx <- dropNA2matrix(assay(tf2peak_cor.se,"cor"))

genes.to.plot <- c("FOXA2","FOXC1","FOXC2")
# i <- "FOXA2"

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    peak = rownames(pvalue.mtx),
    cor = cor.mtx[,i],
    pvalue = pvalue.mtx[,i]
    ) %>% 
    .[!is.na(pvalue)] %>% 
    # .[,qvalue:=p.adjust(pvalue, method = "fdr")] %>%
    # .[,sig:=abs(cor)>0.25 & pvalue<0.01] %>%
    .[,class:="Not correlated"] %>%
    .[cor>0.25 & pvalue<0.50,class:="Positively correlated"] %>%
    .[cor<(-0.25) & pvalue<0.50,class:="Negatively correlated"] %>%
    .[,class:=factor(class, levels=c("Not correlated","Positively correlated","Negatively correlated"))]
  
  p1 <- gghistogram(to.plot, x="cor", y="..density..", fill="gray90", bins=50) +
    labs(x="Correlation coefficient", y="Density") +
    theme(
      axis.text = element_text(size=rel(0.75))
    )
  
  p2 <- ggbarplot(to.plot[,.N,by="class"], x="class", y="N", fill="class") +
    geom_text(aes(x = class, y = N, label = paste("N =",N,"\n"))) +
    scale_x_discrete(drop=F) +
    scale_fill_brewer(palette="Dark2") +
    labs(x="", y="Number of peaks") +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=rel(0.75))
    )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)
  
  pdf(sprintf("%s/%s_histogram_barplot_correlation_events.pdf",io$outdir,i), width = 9, height = 5)
  print(p2)
  dev.off()
  
}

#############
## Heatmap ##
#############

genes.to.plot <- colnames(tf2peak_cor.se)

i <- "GATA1"
for (i in genes.to.plot) {
  
  outfile <- sprintf("%s/heatmaps/%s_heatmap_corr_peaks.pdf",io$outdir,i)
  if (!file.exists(outfile)) {
    print(i)
    
    # peaks.to.plot <- assay(tf2peak_cor.se[,i])[,1][assay(tf2peak_cor.se[,i])[,1]!=0] %>% sort
    peaks.to.plot <- assay(tf2peak_cor.se[,i])[,1][dropNA2matrix(assay(tf2peak_cor.se[,i], "pvalue"))[,1]<0.10] %>% sort
    if (length(peaks.to.plot)>=5) {
      
      to.plot.rna <- logcounts(rna.sce[i,])
      to.plot.acc <- assay(atac.peakMatrix.se[names(peaks.to.plot),])
      
      
      # Filter peaks
      # to.plot.acc <- to.plot.acc[apply(to.plot.acc,1,var)>0.001,]
      # peaks.to.plot <- peaks.to.plot[names(peaks.to.plot)%in%rownames(to.plot.acc)]
      
      # Prepare side information
      annotation_col.df <- data.frame(celltype=opts$celltypes); rownames(annotation_col.df) <- opts$celltypes
      mycolors <- list("celltype"=opts$celltype.colors[opts$celltypes])
      
      annotation_row.df <- data.frame(cor=peaks.to.plot)
      
      to.plot.acc[to.plot.acc>1] <- 1
      
      # Plot heamtap
      p.rna <- pheatmap(
        mat = to.plot.rna, 
        cluster_cols = F, cluster_rows = F,
        color = colorRampPalette(c("gray80", "purple"))(100),
        show_colnames = FALSE, show_rownames = FALSE,
        annotation_col = annotation_col.df,
        annotation_colors = mycolors,
        legend = F, annotation_legend = F,
        silent = TRUE
      )
      
      p.atac <- pheatmap(
        mat = to.plot.acc, 
        cluster_cols = F, cluster_rows = T,
        annotation_row = annotation_row.df,
        show_colnames = FALSE, show_rownames = FALSE,
        legend = F, annotation_legend = F,
        annotation_colors =  list("cor"=colorRampPalette(c("red", "blue"))(10)),
        treeheight_row = 0,
        silent = TRUE
      )
      
      
      p <- cowplot::plot_grid(p.rna$gtable, p.atac$gtable, nrow = 2, rel_heights = c(1/10,9/10))
      pdf(outfile, width=7, height=10)
      print(p)
      dev.off()
      
      # p <- cowplot::plot_grid(p.rna$gtable, p.atac$gtable, nrow = 2, rel_heights = c(1/10,9/10))
      # pdf(sprintf("%s/%s_test.pdf",io$outdir,i), width=5, height=10)
      # print(p)
      # dev.off()
    }
  } else {
    print(sprintf("%s already exists...",outfile))
  }
}
