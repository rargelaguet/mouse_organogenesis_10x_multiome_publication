source(here::here("settings.R"))
source(here::here("utils.R"))

library(GenomicRanges)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(furrr)

#####################
## Define settings ##
#####################

# I/O
# io$peak_annotation <- file.path(io$basedir,"processed/atac/archR/Annotations/CISBP-Scores.rds")
# io$motif2gene <- file.path(io$basedir,"processed/atac/archR/Annotations/CISBP_motif2gene.txt.gz")
# io$feature_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
io$diff_results <- "/Users/rargelaguet/data/fibros_reik_lab/results/atac/archR/deeptools_pseudobulk/ATAC_peaks/differential/25F_vs_63F.tsv"
# io$outdir <- file.path(io$basedir,"results/atac/archR/differential/metacells/alias/PeakMatrix/parsed/pdf/motif_enrichment")
# dir.create(io$outdir, showWarnings = F)


###########################
## Load feature metadata ##
###########################

# feature_metadata.dt <- fread(io$feature_metadata)

############################
## Load motif annotations ##
############################

# peak_annotation.se <- readRDS(io$peak_annotation)

# motif2gene.dt <- fread(io$motif2gene)

###############################
## Load differential results ##
###############################

diff.dt <- fread(io$diff_results) %>%
  .[,sign:=ifelse(diff>0,"Open","Closed")]
table(diff.dt$chr)

# convert to genomicranges
diff.gr <- makeGRangesFromDataFrame(diff.dt, keep.extra.columns = T)

#####################
## Run motifmatchr ##
#####################

# Get PWM List adapted from chromVAR!
data("human_pwms_v2")
motifs <- human_pwms_v2
obj <- .summarizeChromVARMotifs(motifs)
motifs <- obj$motifs
motifSummary <- obj$motifSummary
stopifnot(names(motifs)==rownames(motifSummary))

#########################
## Get motif positions ##
#########################

motifScores <- matchMotifs(
  # pwms = motifs["T_789"],
  pwms = motifs,
  subject = diff.gr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  out = "scores",
  p.cutoff = 1e-3,
  w = 7
)

# tmp <- motifPositions[[1]]
# tmp[seqnames(tmp)=="chr15" & start(tmp)>36706260 & end(tmp)<=36706860]
# tmp[tmp$score>=10]

class(assay(motifScores,"motifScores"))
class(assay(motifScores,"motifMatches"))
class(assay(motifScores,"motifCounts"))

###########################################################
## Run motif enrichment for each differential comparison ##
###########################################################

motifs.to.use <- names(motifs)# %>% head(n=3)

peak_annnotation_match.mtx <- assay(motifScores,"motifMatches")
rownames(peak_annnotation_match.mtx) <- rowData(motifScores)$feature

i <- "25F"; j <- "63F"

# m <- motifs.to.use[1]
motif_enrichment.list <- motifs.to.use %>% map(function(m) { c("Open","Closed") %>% map(function(k) { 
  foreground.peaks <- diff.dt[sign==k & abs(diff)>0.05,feature]
  background.peaks <- diff.dt[abs(diff)<0.01,feature]
  
  foreground.nmatches <- peak_annnotation_match.mtx[foreground.peaks,m][] %>% sum
  background.nmatches <- peak_annnotation_match.mtx[background.peaks,m][] %>% sum
  
  p.value <- phyper(
    foreground.nmatches-1, 
    background.nmatches, 
    length(background.peaks)-background.nmatches,
    length(foreground.peaks), 
    lower.tail = F
  )
  
  data.table(motif=m, groupA=i, groupB=j, sign=k, pval=format(p.value,digits=3))
}) %>% rbindlist })


motif_enrichment.dt <- rbindlist(motif_enrichment.list)
motif_enrichment.dt$pval <- as.numeric(motif_enrichment.dt$pval)

fwrite(motif_enrichment.dt, file.path(io$outdir,"motif_enrichment_diff_acc.txt.gz"), sep="\t", na="NA", quote=F)


##########
## Plot ##
##########

to.plot <- motif_enrichment.dt[sign=="Closed"] %>% setorder(pval) %>% head(n=25) %>%
  .[,tf:=strsplit(motif,"_") %>% map_chr(1)] %>%
  .[,log_pval:=-log10(pval)] %>% .[,tf:=factor(tf, levels=tf)]

ggplot(to.plot, aes_string(x="log_pval", y="tf")) +
  geom_point(aes(size=log_pval)) +
  geom_segment(aes(yend=tf), xend=0, linewidth=0.20) +
  scale_size_continuous(range = c(1.25,3.75)) +
  labs(y="", x="Motif enrichment (p-val)", title="25F vs 63F") +
  theme_classic() +
  guides(size="none") +
  theme(
    plot.title = element_text(size=rel(1.25), hjust=0.5),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(1.20), color="black")
  )
