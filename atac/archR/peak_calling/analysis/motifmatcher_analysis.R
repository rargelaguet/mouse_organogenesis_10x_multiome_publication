here::i_am("atac/archR/peak_calling/analysis/motifmatcher_analysis.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$basedir <- file.path(io$basedir,"test")
io$atac_peak_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
io$motifmatcher_scores <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/CISBP-Scores.rds"))
io$motifmatcher_positions <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/CISBP-Positions.rds"))
io$motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/CISBP_motif2gene.txt.gz"))
io$outdir <- paste0(io$basedir,"/results/atac/archR/peak_calling/motifmatchr"); dir.create(io$outdir, showWarnings = F)

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(io$atac_peak_metadata) %>%
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)]

#####################
## Load motif2gene ##
#####################

motif2gene.dt <- fread(io$motif2gene)

###############################
## Load motifmatcher results ##
###############################

motifmatcher_scores.se <- readRDS(io$motifmatcher_scores)
motifmatcher_positions.se <- readRDS(io$motifmatcher_positions)

################################
## Parse motifmatcher results ##
################################

# Subset motifs
motifs <- intersect(motif2gene.dt$motif, names(motifmatcher_positions.se))
motifmatcher_positions.se <- motifmatcher_positions.se[motifs]
motif2gene.dt <- motif2gene.dt[motif%in%motifs]

# Rename TFs
tmp <- motif2gene.dt$gene; names(tmp) <- motif2gene.dt$motif
names(motifmatcher_positions.se) <- tmp[names(motifmatcher_positions.se)]
colnames(motifmatcher_scores.se) <- tmp[colnames(motifmatcher_scores.se)]
motifmatcher_positions.se <- motifmatcher_positions.se[!duplicated(names(motifmatcher_positions.se))]
motifmatcher_scores.se <- motifmatcher_scores.se[!duplicated(colnames(motifmatcher_scores.se))]

###############################################################
## Plot distribution of motifmatchr scores across all motifs ##
###############################################################

# Create long data.table
tmp <- motifmatcher_positions.se[seqnames(motifmatcher_positions.se)%in%c("chr1")]
to.plot <- unlist(tmp) %>% as.data.table %>% 
  .[,motif:=as.factor(names(unlist(tmp)))] %>% 
  setnames("seqnames","chr") %>% .[,motif_location:=sprintf("%s:%s-%s",chr,start,end)]

p <- ggdensity(to.plot, x="score", fill="strand") +
  labs(x="motifmarchr scores") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.8))
  )
pdf(paste0(io$outdir,"/histogram_motifmatchr_scores_allmotifs.pdf"), width = 6, height = 5)
print(p)
dev.off()

#######################################################
## Plot distribution of motifmatchr scores per motif ##
#######################################################

max.score <- max(motifmatcher_positions.dt$score)
min.score <- min(motifmatcher_positions.dt$score)

for (i in names(motifmatcher_positions.se)) {
  
  print(i)
  to.plot <- motifmatcher_positions.dt[motif==i]
  
  p <- ggdensity(to.plot, x="score", fill="strand") +
    coord_cartesian(xlim=c(min.score-0.1,max.score+0.1)) +
    labs(x=sprintf("%s motifmatchr scores (N=%d)",i,nrow(to.plot))) +
    theme(
      legend.title = element_blank(),
      axis.text = element_text(size=rel(0.8))
    )
  
  pdf(sprintf("%s/histograms_per_TF/%s_histogram_motifmatchr_scores.pdf",io$outdir,i), width = 6, height = 5)
  print(p)
  dev.off()
}


##########################################
## Plot number of binding events per TF ##
##########################################

to.plot <- data.table(
  TF = colnames(motifmatcher_scores.se),
  N = colSums(assay(motifmatcher_scores.se,"motifMatches"))
) %>% setorder(-N) %>% .[,TF:=factor(TF)] %>% 
  .[,log_N:=log10(N)]

p <- ggboxplot(to.plot, x="TF", y="log_N") +
  labs(x="Number of motif matches per TF") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.8))
  )
pdf(paste0(io$outdir,"/boxplot_number_of_motifs_per_TF.pdf"), width = 6, height = 5)
print(p)
dev.off()

p <- gghistogram(to.plot, x="log_N", fill="gray70") +
  labs(x="Number of motif matches per TF (log10)") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.8))
  )
pdf(paste0(io$outdir,"/histogram_number_of_motifs_per_TF.pdf"), width = 6, height = 5)
print(p)
dev.off()

##################################################
## Plot distribution of binding events per peak ##
##################################################

to.plot <- data.table(
  peak = rownames(motifmatcher_scores.se),
  N = rowSums(assay(motifmatcher_scores.se,"motifMatches"))
) %>% merge(peak_metadata.dt[,c("peak","peakType","GC","score")], by=c("peak"))

p <- ggdensity(to.plot, x="N", fill="peakType") +
  labs(x="Number of TF motifs per peak") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.8))
  )

pdf(paste0(io$outdir,"/histogram_number_of_motifs_per_peak.pdf"), width = 6, height = 5)
print(p)
dev.off()


##########
## Save ##
##########

# to.save <- data.table(
#   motif = colnames(motifmatcher_scores.se),
#   N = colSums(assay(motifmatcher_scores.se))
# ) %>% 
#   merge(unique(motifmatcher_positions.dt[,c("width","motif")]),by="motif") %>%
#   setorder(-N)
# fwrite(to.save, paste0(io$outdir,"/number_motifmatches_per_TF.txt.gz"), sep="\t", quote=F)
# 
# to.save <- data.table(
#   peak = rownames(motifmatcher_scores.se),
#   N = rowSums(assay(motifmatcher_scores.se))
# ) %>% merge(peak_metadata.dt[,c("peak","peakType","GC","score")], by=c("peak")) %>%
#   setnames("score","peakScore")
# fwrite(to.save, paste0(io$outdir,"/number_motifmatches_per_peak.txt.gz"), sep="\t", quote=F)

#################################################
## Plot motif score vs number of motif matches ##
#################################################

to.plot <- seq(1,10,by=0.5) %>%  map( function(x) {
    data.table(
      min_score = x,
      peak = rownames(motifmatcher_scores.se),
      N = rowSums(assay(motifmatcher_scores.se,"motifScores")>=x)
    )
  }) %>% rbindlist

to.plot2 <- to.plot[,.(mean=mean(N), sd=sd(N)),by=c("min_score")]

p <- ggplot(to.plot2, aes(x=min_score, y=mean)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.25, colour="black", size=0.5) +
  labs(x="Minimum motif score", y="Total number of motif matches") +
  theme_classic() +
  theme(
    axis.text = element_text(color="black", size=rel(1)),
    axis.title = element_text(color="black", size=rel(1.1))
  )

pdf(paste0(io$outdir,"/motifscore_vs_number_matches_per_peak.pdf"), width = 6, height = 5)
print(p)
dev.off()

##############################
## Inspect individual peaks ##
##############################

number_motifs_per_peak <- rowSums(assay(motifmatcher_scores.se,"motifMatches"))
opts$min_motif_score <- 3
peaks.to.plot <- number_motifs_per_peak[number_motifs_per_peak>=50 & number_motifs_per_peak<=100] %>% names %>% sample(size=25)

for (i in peaks.to.plot) {
  
  peak.gr <- rowRanges(motifmatcher_scores.se[i,])
  
  # TFs <- which(assay(motifmatcher_scores.se[i,],"motifMatches")[1,]) %>% names
  TFs <- which(assay(motifmatcher_scores.se[i,],"motifScores")[1,]>=4) %>% names

  ## START TEST
  # j <- "SREBF1"
  # foo <- motifmatcher_positions.se[[j]][seqnames(motifmatcher_positions.se[[j]])==seqnames(peak.gr)]
  # foo[(seqnames(foo)==seqnames(peak.gr)) & (start(foo)>=start(peak.gr)) & (end(foo)<=end(peak.gr))]
  ## END TEST
  
  motif_locations.dt <- TFs %>% map(function(j) {
    tmp <- motifmatcher_positions.se[[j]][seqnames(motifmatcher_positions.se[[j]])==seqnames(peak.gr)]
    hits <- findOverlaps(
      query = peak.gr,
      subject = tmp,
      ignore.strand = T
    ) %>% subjectHits()
    foo <- tmp[hits]
    foo$TF <- j
    return(foo)
  }) %>% as(., "GRangesList") %>% unlist %>%
    as.data.table %>% setorder(start)
  
  # Plot
  to.plot <- motif_locations.dt %>% .[score<=opts$min_motif_score,score:=opts$min_motif_score]
  # to.plot.text <- to.plot %>% setorder(-score) %>% head(n=15)
  if (nrow(to.plot)>50) {
    to.plot.text <- to.plot[sample(.N,50)]
  } else {
    to.plot.text <- to.plot
  }
  
  
  p <- ggplot(to.plot, aes_string(x="start", y="score")) +
    geom_jitter(aes(size=score, shape=strand), width = 10, height=0.05, fill="gray70") +
    # scale_x_continuous(limits=c(start(peak.gr),end(peak.gr))) +
    # scale_fill_gradientn(colours = terrain.colors(10)) +
    # scale_fill_gradient2(low = "gray50", high = "gray10") +
    scale_shape_manual(values=c(21,24)) +
    scale_size_continuous(range=c(1,4)) +
    ggrepel::geom_text_repel(aes(label=TF), size=3, max.overlaps=Inf, data=to.plot.text) +
    # geom_segment(aes_string(xend="score"), size=0.5, yend=0) +
    labs(x="Genomic position (ATAC peak)", y="Motif score", title=i) +
    theme_classic() +
    guides(fill="none", size="none") +
    theme(
      plot.title = element_text(hjust=0.5),
      legend.position = "none",
      axis.text = element_text(size=rel(0.8), color="black"),
      axis.title = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/peak_%s_motifs_location.pdf",io$outdir,sub(":","-",i)), width=7, height=5)
  print(p)
  dev.off()
}


