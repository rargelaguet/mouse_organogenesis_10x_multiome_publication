
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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/per_peak"); dir.create(io$outdir, showWarnings = F)

# Options
opts$motif_annotation <- "Motif_cisbp_lenient"

###############
## Load data ##
###############

# io$tf2peak_cor.se <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
io$tf2peak_cor.se <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment_lenient.rds")
tf2peak_cor.se <- readRDS(io$tf2peak_cor.se)

# io$tf2peak_cor.dt <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_lenient.txt.gz")
# cor_dt <- fread(io$tf2peak_cor.dt) %>%
#   .[!is.na(cor)] %>%
#   .[,cor_sign:=c("-","+")[(cor>0)+1]]

######################
## Load motifmatchR ##
######################

io$motifmatcher.se <- sprintf("%s/Annotations/%s-Matches-In-Peaks.rds",io$archR.directory,opts$motif_annotation)
io$motifmatcher_positions.se <- sprintf("%s/Annotations/%s-Positions-In-Peaks.rds",io$archR.directory,opts$motif_annotation)

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
  source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR_positions.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
  source("/homes/ricard/gastrulation_multiome_10x/load_motifmatchR_positions.R")
} else {
  stop("Computer not recognised")
}

#####################################################################
## Plot histogram of the nº of TF motifs per peak before and after ##
#####################################################################

# (TO-DO: USE MINIMUM SCORE THRESHOLD)
opts$min.corr <- 0.50

before_cor_filtering.dt <- data.table(
  peak = rownames(motifmatcher.se),
  N = rowSums(assay(motifmatcher.se))
)

after_cor_filtering.dt <- data.table(
  peak = rownames(tf2peak_cor.se),
  # N = rowSums(dropNA2matrix(assay(tf2peak_cor.se,"pvalue"))<0.01,na.rm=T)
  # N = rowSums(abs(dropNA2matrix(assay(tf2peak_cor.se,"cor")))>=opts$min.corr,na.rm=T)
  N = Matrix::rowSums(abs(assay(tf2peak_cor.se,"cor"))>=opts$min.corr,na.rm=T)
) 

to.plot <-  merge(before_cor_filtering.dt, after_cor_filtering.dt, by="peak", suffixes=c("_chromVAR","_chromRNA")) %>%
  melt(id.vars=c("peak"))

p <- ggdensity(to.plot, x="value", fill="variable", y="..density..") +
  geom_vline(aes(xintercept=V1), color="black", data=to.plot[,median(value),by="variable"], linetype="dashed") +
  labs(x="Number of TF motifs per peak", y="Density") +
  scale_x_continuous(limits = c(0,165)) +
  scale_fill_brewer(palette="Dark2") +
  # guides(color=F) + 
  theme(
    axis.title = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.text.x = element_text(size=rel(0.80), color="black"),
    axis.text.y = element_text(size=rel(0.50), color="black")
  )


pdf(sprintf("%s/density_number_motifs_per_peak_before_and_after.pdf",io$outdir), width=5, height=4)
print(p)
dev.off()
  
###########################
## Plot individual peaks ##
###########################

number_correlated_TFs_per_peak <- rowSums(abs(assay(tf2peak_cor.se))>0.5)

# number_motifs_per_peak <- rowSums(assay(tf2peak_cor.se))
opts$min.score <- 4

set.seed(42)
peaks.to.plot <- number_correlated_TFs_per_peak[number_correlated_TFs_per_peak>=10 & number_correlated_TFs_per_peak<=50] %>% names %>% sample(size=25)
# peaks.to.plot <- "chr1:3035578-3036178"

for (i in peaks.to.plot) {
  
  # Identify motif locations + scores
  peak.gr <- rowRanges(motifmatcher.se[i,])
  TFs <- assay(motifmatcher.se[i])[1,] %>% which %>% names
  motif_locations.dt <- TFs %>% map(function(j) {
    hits <- findOverlaps(
      query = peak.gr,
      subject = motifmatcher_positions.se[[j]]
    ) %>% subjectHits()
    tmp <- motifmatcher_positions.se[[j]][hits]
    tmp$TF <- j
    return(tmp)
  }) %>% as(., "GRangesList") %>% unlist %>%
    as.data.table %>% setorder(start)
  
  # Add correlation with TFexpression
  tf2peak_cor_tmp.dt <- data.table(
    TF = colnames(tf2peak_cor.se),
    cor = assay(tf2peak_cor.se[i,])[1,]
  )
  to.plot <- motif_locations.dt %>% merge(tf2peak_cor_tmp.dt,by="TF") %>% .[,abs_cor:=abs(cor)]
  
  # Remove palindromic motif duplicates
  to.plot[,N:=.N,by=c("seqnames","start","end","TF")] 
  to.plot <- rbind(to.plot[N==1],to.plot[N==2 & strand=="+"]) %>% .[,N:=NULL]
  
  # Filter by mininum score
  to.plot <- to.plot %>% .[score>=opts$min.score]
  # to.plot <- motif_locations.dt %>% .[score<=opts$min.score,score:=opts$min.score]
  
  # Plot
  # to.plot.text <- to.plot %>% setorder(-score) %>% head(n=15)
  if (nrow(to.plot)>25) {
    to.plot.text <- to.plot[score>7] %>% setorder(abs_cor) %>% tail(n=25)
  } else {
    to.plot.text <- to.plot
  }
  
  p <- ggplot(to.plot, aes_string(x="start", y="score")) +
    geom_jitter(aes(fill=cor, size=score, shape=strand, alpha=abs_cor), width = 5, height=0.1) +
    scale_x_continuous(limits=c(start(peak.gr),end(peak.gr))) +
    # scale_fill_gradientn(colours = rev(terrain.colors(10))) +
    # scale_fill_gradient2(low = "gray50", high = "gray10") +
    scale_fill_gradient2(low = "blue", mid="gray90", high = "red", limits=c(-1,1)) +
    scale_shape_manual(values=c(21,24)) +
    scale_size_continuous(range=c(0.5,3)) +
    scale_alpha_continuous(range=c(0.25,1)) +
    ggrepel::geom_text_repel(aes(label=TF), size=3, max.overlaps=Inf, data=to.plot.text) +
    # geom_segment(aes_string(xend="score"), size=0.5, yend=0) +
    labs(x="Genomic position (ATAC peak)", y="Motif score", title=i) +
    theme_classic() +
    guides(size=F, alpha=F) +
    theme(
      plot.title = element_text(hjust=0.5),
      legend.position = "top",
      axis.text = element_text(size=rel(0.8), color="black"),
      axis.title = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/peak_%s_rna_vs_acc_TFmotifs_location.pdf",io$outdir,sub(":","-",i)), width=6, height=4)
  print(p)
  dev.off()
}

# Plot legend
# p <- ggplot(to.plot, aes_string(x="start", y="score")) +
#   geom_point(aes(shape=strand), size=5, color="gray10") +
#   # scale_shape_manual(values=c(1,2)) +
#   # guides(fill=F, size=F) +
#   theme_classic() +
#   theme(
#     legend.position = "top"
#   )
# 
# pdf(sprintf("%s/legend.pdf",io$outdir), width=7, height=5)
# print(p)
# dev.off()

###############################################################
## Co-occurence analysis ##
###############################################################



#############################################################
## IGNORE BELOW ##
#############################################################


to.plot <- cor_dt %>%
  merge(peak_metadata.dt[,c("peak","peakType")]) %>%
  .[,.N,by=c("peakType","TF")] 

TFs.to.plot <- cor_dt[,.N,by="TF"] %>% setorder(-N) %>% head(n=100) %>% .$TF
to.plot <- to.plot[TF%in%TFs.to.plot] %>% .[,TF:=factor(TF,levels=TFs.to.plot)]

p <- ggbarplot(to.plot, x="TF", y="N", fill="peakType") +
  coord_flip() +
  # scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of correlation events with peaks") +
  theme(
    
    axis.text.x = element_text(colour="black",size=rel(0.8)),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )


pdf(sprintf("%s/barplot_number_correlation_events_perTF_perPeakType.pdf",io$outdir), width = 6, height = 12)
print(p)
dev.off()

