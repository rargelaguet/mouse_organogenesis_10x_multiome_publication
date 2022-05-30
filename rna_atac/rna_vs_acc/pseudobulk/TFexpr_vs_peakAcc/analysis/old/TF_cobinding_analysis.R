
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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cobinding")

calculate_cobinding_stats <- function(mtx) {
  pval.mtx <- matrix(NA,nrow = nrow(mtx), ncol=ncol(mtx))
  dimnames(pval.mtx) <- dimnames(mtx)
  for (i in 1:nrow(mtx)) {
    m <- sum(mtx[i,])
    for (j in i:ncol(mtx)) {
      if (i!=j) {
        n <- sum(mtx[j,])
        expected = (m*n)/N
        observed = mtx[i,j]
        pval.mtx[i,j] <- pval.mtx[j,i] <-poisson.test(observed, r = expected, alternative = c("greater"))$p.value
      }
    }
  }
  return(pval.mtx)
}

###############
## Load data ##
###############

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
tf2peak_cor.se <- readRDS(io$file)

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc.txt.gz")
tf2peak_cor.dt <- fread(io$file) %>%
  .[!is.na(cor)] %>%
  .[,cor_sign:=c("-","+")[(cor>0)+1]]


##################################################
## Load pseudobulk ATAC and RNA expression data ##
##################################################

opts$motif_annotation <- "Motif_cisbp"
io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}


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

#######################
## Load marker peaks ##
#######################

marker_peaks.dt <- fread(io$markers_peaks)

#####################################
## Create motif2peak binary matrix ##
#####################################

# Create binding matrix
# motif2peak.mtx <- cor_dt %>%
#   dcast(peak~TF, fun.aggregate=function(x){as.integer(length(x) > 0)}) %>% 
#   matrix.please %>%
#   as("lgCMatrix")


#############################
## Detect cobinding events ##
#############################

# Subset to peaks with N>=2
motif2peak.mtx.filt <- motif2peak.mtx[rowSums(motif2peak.mtx)>=2,] %>% as.matrix

foo <- crossprod(motif2peak.mtx.filt)
foo <- foo[rowSums(foo)>=5,colSums(foo)>=5]
bar <- sweep(foo,1,apply(foo,1,max),FUN = "/")

################
## Statistics ##
################

# i <- "EOMES"
# j <- "T"
# m <- sum(foo[i,])
# n <- sum(foo[j,])
# N <- nrow(motif2peak.mtx.filt)
# expected = (m*n)/N
# observed = foo[i,j]
# poisson.test(observed, T = 1, r = expected, alternative = c("greater"))$p.value

pval.mtx <- matrix(NA,nrow = nrow(foo), ncol=ncol(foo))
dimnames(pval.mtx) <- dimnames(foo)
for (i in 1:nrow(foo)) {
  m <- sum(foo[i,])
  for (j in 1:ncol(foo)) {
    n <- sum(foo[j,])
    expected = (m*n)/N
    observed = foo[i,j]
    pval.mtx[i,j] <- poisson.test(observed, r = expected, alternative = c("greater"))$p.value
  }
}
diag(pval.mtx) <- NA

###################################
## Plot cobinding results per TF ##
###################################

cobinding.dt <- pval.mtx %>% as.data.table(keep.rownames = T) %>%
  melt(id.vars="rn") %>% setnames(c("TF1","TF2","pvalue")) %>%
  .[!is.na(pvalue)] %>%
  .[,log_padj_fdr := -log(p.adjust(pvalue+1e-50, method="fdr")), by="TF1"]

tfs.to.plot <- cobinding.dt %>% .[TF1!=TF2,max(log_padj_fdr),by="TF1"] %>% .[V1>2,TF1]

for (i in tfs.to.plot) {
  
  to.plot <- cobinding.dt[TF1==i] %>% 
    setorder(-log_padj_fdr) %>% head(n=15) %>%
    .[,TF2:=factor(TF2,levels=TF2)]
  
  p <- ggplot(to.plot, aes_string(x="TF2", y="log_padj_fdr"), fill="gray70") +
    geom_point(size=2) +
    geom_segment(aes_string(xend="TF2"), size=0.5, yend=0) +
    coord_flip() +
    labs(x="", y="Cobinding score (log10 p-value)") +
    theme_classic() +
    theme(
      axis.ticks.y = element_blank(),
      axis.text = element_text(size=rel(0.75), color="black")
    )
  
  pdf(sprintf("%s/cobinding/%s_cobinding.pdf",io$outdir,i), width = 5, height = 4)
  print(p)
  dev.off()
}


##############################################
## Detect cobinding events for each lineage ##
##############################################

# i <- "Epiblast"
for (i in unique(marker_peaks.dt$celltype)) {
  
  dir.create(sprintf("%s/%s",io$outdir,i), showWarnings = F)
  motif2peak.mtx.filt <- motif2peak.mtx[rownames(motif2peak.mtx)%in%marker_peaks.dt[celltype==i,idx],] %>% as.matrix
  
  tmp <- crossprod(motif2peak.mtx.filt)
  tmp <- tmp[rowSums(tmp)>=5,colSums(tmp)>=5]
  
  pval.mtx <- calculate_cobinding_stats(tmp)
  
  # to.plot <- sweep(tmp,1,apply(tmp,1,max),FUN = "/")
  # pheatmap::pheatmap(to.plot, cluster_rows = F, cluster_cols = F)
  
  cobinding.dt <- pval.mtx %>% as.data.table(keep.rownames = T) %>%
    melt(id.vars="rn") %>% setnames(c("TF1","TF2","pvalue")) %>%
    .[!is.na(pvalue)] %>%
    .[,log_padj_fdr := -log(p.adjust(pvalue+1e-50, method="fdr")), by="TF1"]
  
  tfs.to.plot <- cobinding.dt %>% .[TF1!=TF2,max(log_padj_fdr),by="TF1"] %>% .[V1>5,TF1]
  
  for (j in tfs.to.plot) {
    
    to.plot <- cobinding.dt[TF1==j] %>% 
      setorder(-log_padj_fdr) %>% head(n=15) %>%
      .[,TF2:=factor(TF2,levels=TF2)]
    
    p <- ggplot(to.plot, aes_string(x="TF2", y="log_padj_fdr"), fill="gray70") +
      geom_point(size=2) +
      geom_segment(aes_string(xend="TF2"), size=0.5, yend=0) +
      coord_flip() +
      labs(x="", y="Cobinding score (log10 p-value)") +
      theme_classic() +
      theme(
        axis.ticks.y = element_blank(),
        axis.text = element_text(size=rel(0.75), color="black")
      )
    
    pdf(sprintf("%s/%s/%s_%s_cobinding.pdf",io$outdir,i,i,j), width = 5, height = 4)
    print(p)
    dev.off()
  }
  
}


#############
## Explore ##
#############



asd <- tf2peak_cor.dt %>% copy %>%
  .[,Ntotal:=.N,by="TF"] %>%
  .[Ntotal>15] %>%
  .[,.N,by = c("TF","cor_sign","Ntotal")] %>%
  dcast(TF~cor_sign,value.var="N", fill = 0) %>%
  .[,ratio:=(`+`+0.1)/(`-`+0.1)] %>%
  .[,class:=c("Repressor","Activator")[(ratio>0.5)+1]] %>%
  .[ratio>0.25 & ratio<0.75,class:="Mixed"]


i <- "SOX4"
atac.peakMatrix.filt.se <- atac.peakMatrix.se[peaks,]

all.peaks <- which(!is.na(dropNA2matrix(assay(tf2peak_cor.se[,i],"cor"))[,1])) %>% names
corr_pos.peaks <- which(dropNA2matrix(assay(tf2peak_cor.se[,i],"cor"))[,1]>=0.50) %>% names
corr_neg.peaks <- which(dropNA2matrix(assay(tf2peak_cor.se[,i],"cor"))[,1]<=(-0.50)) %>% names
corr.peaks <- which(dropNA2matrix(assay(tf2peak_cor.se[,i],"cor"))[,1]>=0.50) %>% names
noncorr.peaks <- all.peaks[!all.peaks%in%corr.peaks]

motifmatcher_corrpeaks.se <- motifmatcher.se[corr.peaks,]
motifmatcher_corrnegpeaks.se <- motifmatcher.se[corr_neg.peaks,]
motifmatcher_corrpospeaks.se <- motifmatcher.se[corr_pos.peaks,]
motifmatcher_noncorrpeaks.se <- motifmatcher.se[noncorr.peaks,]

foo <- Matrix::colMeans(assay(motifmatcher_corrnegpeaks.se))
bar <- Matrix::colMeans(assay(motifmatcher_corrpospeaks.se))

# to.plot <- data.table(
#   TF = colnames(motifmatcher.se),
#   fraction_corr_peaks = foo,
#   fraction_noncorr_peaks =   bar
# ) %>% .[,diff:=fraction_corr_peaks-fraction_noncorr_peaks] %>% setorder(-diff)

to.plot <- data.table(
  TF = colnames(motifmatcher.se),
  fraction_pos_peaks = foo,
  fraction_neg_peaks =   bar
) %>% .[,diff:=fraction_pos_peaks-fraction_neg_peaks] %>% setorder(-diff) %>%
  .[,abs_diff:=abs(diff)] %>% .[abs_diff>=0.075]

to.plot2 <- to.plot %>% merge(asd, by="TF")
to.plot2 %>% setorder(-diff) %>% .[,TF:=factor(TF,levels=TF)]


p <- ggplot(to.plot2, aes(x=TF, y=diff, color=class)) +
  geom_segment(aes_string(xend="TF"), size=0.5, yend=0) +
  geom_point() +
  scale_x_discrete(drop=FALSE) +
  coord_flip() +
  labs(y="Difference in the proportion of motifs") +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(1.2)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(1.1), color="black"),
    axis.text.x = element_text(size=rel(1.1), color="black")
  )

# pdf(sprintf("%s/per_class/barplots_%s.pdf",io$outdir,i), width=7, height=7)
print(p)
# dev.off()
