library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$basedir <- file.path(io$basedir,"test")
# io$pwm <- file.path(io$basedir,"results/atac/archR/motif_seqlogo/PWMatrixList.rds")
io$chip_dir.prefix <- "/Users/argelagr/data/mm10_regulation/TF_ChIP"
io$virtual_chip_dir <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/pseudobulk/CISBP")
io$motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s_motif2gene.txt.gz",opts$motif_annotation))
io$chip_peaks.files <- c(
  "CDX2" = file.path(io$chip_dir.prefix, "Cdx2_NMP_invitro/fastq/Cdx2_wtEpiSCs_24h_macs2_broad/Cdx2_wtEpiSCs_24h_peaks.broadPeak"),
  "TAL1" = file.path(io$chip_dir.prefix, "TAL1_blood_progenitors/fastq/HP_Tal1_macs2_broad/HP_Tal1_peaks.broadPeak"),
  "GATA1" = file.path(io$chip_dir.prefix, "GATA1_blood_progenitors/original/old/macs2/Gata1_chip_peaks.narrowPeak.gz"),
  # "GATA1" = file.path(io$chip_dir.prefix, "GATA1_blood_progenitors/fastq/HP_Gata1_macs2_broad/HP_Gata1_peaks.broadPeak"),
  # "FOXA2" = file.path(io$chip_dir.prefix, "FOXA2_ESC_endoderm_differentiation/macs2/foxa2_chip_peaks.narrowPeak.gz"),
  "FOXA2" = file.path(io$chip_dir.prefix, "FOXA2_ESC_endoderm_differentiation/original/old/macs2/foxa2_chip_peaks.narrowPeak.gz"),
  "GATA4" = file.path(io$chip_dir.prefix, "GATA4_ESC_mesoderm_differentiation/fastq/Gata4_WT_CP_macs2_broad/Gata4_WT_CP_peaks.broadPeak"),
  "TBX5" = file.path(io$chip_dir.prefix, "TBX5_cardiomyocytes/fastq/Tbx5_WT_CP_macs2_broad/Tbx5_WT_CP_peaks.broadPeak"),
  "RUNX1" = file.path(io$chip_dir.prefix, "RUNX1_blood_progenitors/fastq/HP_Runx1_macs2_broad/HP_Runx1_peaks.broadPeak"),
  "NKX2-5" = file.path(io$chip_dir.prefix, "NKX2-5_cardiomyocytes/fastq/Nkx2-5_WT_CP_macs2_broad/Nkx2-5_WT_CP_peaks.broadPeak")
)
io$outdir <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/pseudobulk/CISBP/validation"); dir.create(io$outdir, showWarnings = F)


# Options
opts$motif_annotation <- "CISBP"
opts$TFs <- names(io$chip_peaks.files)
# opts$TFs <- "FOXA2"

###############
## Load data ##
###############

# Load PWMs
io$peak_annotation_file <- file.path(io$archR.directory,"Annotations/peakAnnotation.rds")
pwms <- readRDS(io$peak_annotation_file)[[opts$motif_annotation]][["motifs"]]

# Load ChIP-seq peak calling results
# (TO-DO: EXTEND PEAKS???)
chip.dt <- opts$TFs %>% map(function(i) {
 fread(io$chip_peaks.files[[i]], select=c(1,2,3,5)) %>%
  setnames(c("chr","start","end","score")) %>%
  .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
  .[!chr%in%c("chrM","chrMT")] %>%
  .[,id:=sprintf("%s:%s-%s",chr,start,end)] %>%
  # .[score>=opts$min_chip_score] %>% 
  setorder(-score) %>% head(n=1000) %>% # We consider the top N ChIP-seq peaks for each sample
  .[,tf:=i] %>%
    return
}) %>% rbindlist %>% setkey(chr,start,end)

chip.dt[,.N,by="tf"]

# Load virtual ChIP-seq
chip_insilico.dt <- opts$TFs %>% map(function(i) {
  fread(sprintf("%s/%s.bed.gz",io$virtual_chip_dir,i)) %>%
    setnames(c("chr","start","end","score")) %>%
    .[,id:=sprintf("%s:%s-%s",chr,start,end)] %>%
    .[score>=0] %>%
    .[,tf:=i] %>%
    return
}) %>% rbindlist %>% setkey(chr,start,end)

# Load motif annotation
motif2gene.dt <- fread(io$motif2gene)

##################
## Prepare data ##
##################

# Rename TFs
# names(pwms)
motif2gene_filt.dt <- motif2gene.dt[gene%in%opts$TFs]
tmp <- motif2gene_filt.dt$gene; names(tmp) <- motif2gene_filt.dt$motif
pwms <- pwms[names(pwms)%in%motif2gene_filt.dt$motif]
names(pwms) <- tmp[names(pwms)]

###################################################
## Find ChIP-seq peaks that contain the TF motif ##
###################################################

opts$peak_extension <- 100

# i <- "TAL1"
chip_peaks_with_motif.dt <- opts$TFs %>% map(function(i) {
  
  gr <- makeGRangesFromDataFrame(chip.dt[tf==i], keep.extra.columns = T)
  
  motifmatcher.se <- motifmatchr::matchMotifs(
    pwms = pwms[i],
    subject = gr+opts$peak_extension,  
    genome = BSgenome.Mmusculus.UCSC.mm10, 
    out = "matches", 
    p.cutoff = 0.001, 
    w = 7
  ); rownames(motifmatcher.se) <- gr$id
  
  peaks.with.motif <- names(which(assay(motifmatcher.se)[,1]))
  return(chip.dt[tf==i & id%in%peaks.with.motif])
  
}) %>% rbindlist %>% setkey(chr,start,end)

###############################
## Plot TF motif specificity ##
###############################

foo <- chip_peaks_with_motif.dt[,.N,by="tf"] %>% setnames("N","peaks_with_motif")
bar <- chip.dt[,.N,by="tf"] %>% setnames("N","total_peaks")

to.plot <- merge(foo,bar,by="tf") %>% 
  .[,percentage_peaks_motif:=peaks_with_motif/total_peaks]# %>%
  # melt(id.vars="tf",value.name="N")

ggbarplot(to.plot, x="tf", y="percentage_peaks_motif", fill="gray70") +
  geom_hline(yintercept=1, linetype="dashed") +
  labs(x="", y="% of peaks with TF motif") +
  theme(
    axis.text.y = element_text(size=rel(0.75))
  )

##############
## Plot ROC ##
#############

# Define range of binding scores
# seq.ranges <- seq(0,max(chip_insilico.dt$score)-0.05,by=0.01)
seq.ranges <- seq(0,0.50,by=0.01)
names(seq.ranges) <- as.character(1:length(seq.ranges))

# ov.dt <- foverlaps(chip_insilico.dt, chip.dt) %>% setnames("i.score","insilico_score")

roc.dt <- opts$TFs %>% map(function(i) {
  
  ov.dt <- foverlaps(chip_insilico.dt[tf==i], chip_peaks_with_motif.dt[tf==i]) %>% setnames("i.score","insilico_score")
  
  seq.ranges %>% map(function(j) {
    
    true_positives = sum(!is.na(ov.dt[insilico_score>=j,start]))
    false_negatives = sum(!is.na(ov.dt[insilico_score<=j,start]))
    false_positives = sum(is.na(ov.dt[insilico_score>=j,start]))
    true_negatives = sum(is.na(ov.dt[insilico_score<=j,start]))
    
    true_positive_rate = true_positives / (true_positives + false_negatives)
    false_positive_rate = false_positives / (false_positives+true_negatives)
    
    data.table(tf=i, min_score=j, tpr=true_positive_rate, fpr=false_positive_rate)
    
  }) %>% rbindlist
}) %>% rbindlist

p <- ggplot(roc.dt, aes_string(x="fpr", y="tpr", color="tf")) +
  scale_color_brewer(palette="Dark2") +
  geom_line(size=1) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(x="False positive rate", y="True positive rate") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(color="black")
  )

pdf(file.path(io$outdir,"virtual_chip_roc_cuves.pdf"), width=8, height=5)
print(p)
dev.off()

