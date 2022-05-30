
# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

library(GenomicRanges)
library(rtracklayer)

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

# I/O
io$basedir <- file.path(io$basedir,"test")
io$chip_dir.prefix <- "~/data/mm10_regulation/TF_ChIP"
io$peak_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
io$virtual_chip.dir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s",opts$motif_annotation))
io$virtual_chip.mtx <- file.path(io$virtual_chip.dir,"virtual_chip_matrix.rds")
io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s/validation",opts$motif_annotation)); dir.create(io$outdir, showWarnings = F)

# TO-DO: TRY NEW BIGWIG FILES
io$chip.files <- c(
  "CDX2" = file.path(io$chip_dir.prefix, "Cdx2_NMP_invitro/bigwig/GSM2253707_Cdx2_wtEpiSCs_24h_rep1_mm10.bw"),
  "TAL1" = file.path(io$chip_dir.prefix, "TAL1_blood_progenitors/bigwig/GSM1692858_HP_Tal1.bw"),
  "GATA1" = file.path(io$chip_dir.prefix, "GATA1_blood_progenitors/original/GSM1692851_HP_Gata1.bw"),
  "FOXA2" = file.path(io$chip_dir.prefix, "FOXA2_ESC_endoderm_differentiation/original/Foxa2.d5FS.bw"),
  "GATA4" = file.path(io$chip_dir.prefix, "GATA4_ESC_mesoderm_differentiation/original/GSM3223330_Gata4.d5FS.r1_m1.ucsc.bigWig"),
  "TBX5" = file.path(io$chip_dir.prefix, "TBX5_cardiomyocytes/fastq/Tbx5_WT_CP.bw")
  
  # TO_IGNORE
  # "RUNX1" = file.path(io$chip_dir.prefix, "RUNX1_blood_progenitors/bigwig/GSM1692856_HP_Runx1.bw"),
  # "NKX2-5" = file.path(io$chip_dir.prefix, "NKX2-5_cardiomyocytes/bigwig/GSM2054327_Nkx_WT_CP_small.bw")
  # "ZIC2" = file.path(io$chip_dir.prefix, "Matsuda2017_POU5F1_ZIC2/ZIC2/GSM1924744_ZIC2peaks_small.bw"),
  # "T" = file.path(io$chip_dir.prefix, "Koch2017_Sox2_T/GSM2454138_ChIP_T_small.bw"),
  # "POU5F1" = "/Users/ricard/Downloads/Matsuda2017/POU5F1/GSM1924747_POU5F1_small.bw",
  # "SOX2" = "/Users/ricard/data/mm10_regulation/TF_ChIP/Koch2017_Sox2_T/GSM2454137_ChIP_Sox2_small.bw",
)


######################
## Load motifmatchR ##
######################

# io$motifmatcher.se <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,opts$motif_annotation)
# motifmatcher.se <- readRDS(io$motifmatcher.se)
# source(here::here("load_motifmatchR.R"))

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(io$peak_metadata) %>%
  .[,c("chr","start","end","score")] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

peak_metadata.gr <- makeGRangesFromDataFrame(peak_metadata.dt, keep.extra.columns = F)

########################
## Load ChIP-seq data ##
########################

stopifnot(file.exists(io$chip.files))
TFs <- names(io$chip.files)

# opts$min.value <- 0.1
# opts$max.value <- 100

# From bigwig files
# chip_dt.list <- TFs %>%
#   map(function(i) {
#     import.bw(BigWigFile(io$chip.files[[i]]), selection = peak_metadata.gr) %>%
#     # import.bw(BigWigFile(io$chip.files[[i]])) %>%
#       as.data.table %>%
#       .[,c(1,2,3,6)] %>% setnames(c("chr", "start", "end","value")) %>%
#       # .[value>=opts$min.value] %>%
#       .[value>=opts$max.value,value:=opts$max.value] %>%
#       # .[, chip_anno := as.factor(i)] %>%
#       setkey(chr,start,end)
#   })# %>% rbindlist %>% setkey(chr,start,end)
# names(chip_dt.list) <- TFs


###################################
## Load virtual ChIP-seq library ##
###################################

# Load matrix
virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)[,TFs]

# temporary fix
# rownames(virtual_chip.mtx) <- sub("-",":",rownames(virtual_chip.mtx))

######################################
## calculate background ChIP signal ##
######################################

# chip_background.dt <- TFs %>% map(function(i) {
#   data.table(
#     chip_anno = i,
#     median = median(import.bw(BigWigFile(io$chip.files[[i]]))$score)
#   )
# }) %>% rbindlist
# fwrite(chip_background.dt, file.path(io$outdir, "chip_background.txt.gz"))

chip_background.dt <- fread(file.path(io$outdir, "chip_background.txt.gz"))
stopifnot(TFs%in%chip_background.dt$chip_anno)

######################
## calculate signal ##
######################

# annoying stuff
peak_metadata_ucsc.gr <- peak_metadata.gr; seqlevelsStyle(peak_metadata_ucsc.gr) <- 'UCSC'
peak_metadata_ncbi.gr <- peak_metadata.gr; seqlevelsStyle(peak_metadata_ncbi.gr) <- 'NCBI'
tf2seqstyle <- rep("UCSC",length(TFs)); names(tf2seqstyle) <- TFs
tf2seqstyle["TBX5"] <- "NCBI"

# i <- "CDX2"
to.plot <- TFs %>% map(function(i) {
  print(i)
  
  # Load virtual ChIP-seq library
  virtual_chip.dt <- fread(sprintf("%s/%s.txt.gz",io$virtual_chip.dir,i)) %>%
    setnames("peak","idx") %>%
    .[,idx:=str_replace(idx,":","-")] %>%
    .[,chr:=strsplit(idx,"-") %>% map_chr(1)] %>%
    .[,start:=as.integer(strsplit(idx,"-") %>% map_chr(2))] %>%
    .[,end:=as.integer(strsplit(idx,"-") %>% map_chr(3))] %>%
    .[,idx:=str_replace(idx,"-",":")] %>% 
    setkey(chr,start,end)
  
  # Load ground truth ChIP-seq data
  if (tf2seqstyle[[i]]=="UCSC") {
    tmp <- peak_metadata_ucsc.gr
  } else if (tf2seqstyle[[i]]=="NCBI") {
    tmp <- peak_metadata_ncbi.gr
  }
  chip.dt <- import.bw(BigWigFile(io$chip.files[[i]]), selection = tmp) %>% as.data.table %>%
    .[,c(1,2,3,6)] %>% setnames(c("chr", "start", "end","value")) %>%
    .[,chr:=as.character(chr)] %>%
    .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
    # .[value>=opts$min.value] %>% .[value>=opts$max.value,value:=opts$max.value] %>% 
  setkey(chr,start,end)
  
  # overlap ATAC peaks with ChIP-seq signal and quantify ChIP signal
  virtual_chip.dt[,c("idx","chr","start","end")] %>%
    foverlaps(chip.dt, nomatch=0) %>% 
    .[,.(value=sum(value)), b=c("idx")] %>%
    merge(virtual_chip.dt[,c("idx","score","correlation_score","max_accessibility_score","motif_score")]) %>%
    .[,chip_anno:=as.factor(i)] %>%
    return
  
}) %>% rbindlist

# fwrite(to.plot, file.path(io$outdir, "to_plot.txt.gz"))
# to.plot <- fread(file.path(io$outdir, "to_plot.txt.gz"))
# stopifnot(TFs%in%to.plot$chip_anno)


######################################################################
## Scatterplots of in silico chip-seq score vs real chip-seq values ##
######################################################################

# max.score <- 0.75
seq.ranges <- seq(0,1,by=0.10); names(seq.ranges) <- as.character(1:length(seq.ranges))

foo <- to.plot %>%
  .[score>=0] %>% 
  merge(chip_background.dt,by="chip_anno") %>%
  .[,log_value:=log(value+1)] %>%
  .[,c("chip_anno","idx","log_value","correlation_score","max_accessibility_score","motif_score","score")] %>%
  setnames("score","score_rna_atac")


to.plot.compare_models <- foo %>%
  # Model 1: DNA
  # .[,score_dna:=minmax.normalisation(motif_score), by="chip_anno"] %>% 
  # Model 2: DNA + ATAC
  .[,score_dna_atac:=minmax.normalisation(max_accessibility_score * motif_score), by="chip_anno"] %>% .[,score_dna_atac:=score_dna_atac+0.001] %>%
  # Model 3: DNA + ATAC + RNA
  .[,score_rna_atac:=minmax.normalisation(score_rna_atac), by="chip_anno"] %>% .[,score_rna_atac:=score_rna_atac+0.001] %>%
  # Discretise values
  # .[,score_dna_group:=cut(abs(score_dna), breaks=seq.ranges)] %>% .[,score_dna_group:=seq.ranges[as.numeric(score_dna_group)]] %>%
  .[,score_dna_atac_group:=cut(abs(score_dna_atac), breaks=seq.ranges)] %>% .[,score_dna_atac_group:=seq.ranges[as.numeric(score_dna_atac_group)]] %>%
  .[,score_rna_atac_group:=cut(abs(score_rna_atac), breaks=seq.ranges)] %>% .[,score_rna_atac_group:=seq.ranges[as.numeric(score_rna_atac_group)]] %>%
  # Prepare for plotting
  melt(id.vars=c("chip_anno","idx","log_value"), measure.vars=c("score_dna_atac_group","score_rna_atac_group"), variable.name="model", value.name="predicted_value") %>%
  .[,N:=.N,by=c("model","chip_anno","predicted_value")] %>%
  .[,.(mean=mean(log_value,na.rm=T), sd=sd(log_value,na.rm=T), N=.N), by=c("model","chip_anno","predicted_value")] %>%
  # Filter settings with a small number of observations
  .[!is.na(mean) & N>=25]


p <- ggline(to.plot.compare_models, x="predicted_value", y="mean", color="model", plot_type="b") +
  facet_wrap(~chip_anno, scales="free_y") +
  # geom_errorbar(aes(ymin=log_value-sd, ymax=log_value-sd), width=.2,) +
  scale_color_brewer(palette="Dark2", labels = c("ATAC", "RNA + ATAC")) +
  # scale_color_discrete(labels = c("DNA", "DNA + ATAC", "DNA + ATAC + RNA")) +
  labs(y="Observed ChIP-seq signal", x="In silico binding score") +
  theme(
    strip.background = element_rect(colour="black", fill=NA),
    axis.text = element_text(size=rel(0.55)),
    axis.title = element_text(size=rel(0.8)),
    legend.title = element_blank(),
    # legend.position = c(.26,.87)
    # legend.position = c(.06,.94)
    legend.position = "top"
  )

pdf(file.path(io$outdir,"insilico_vs_chipseq_lineplot.pdf"), width=6, height=4)
print(p)
dev.off()

##############################################
## Compare positive vs negative predictions ##
##############################################

seq.ranges <- seq(0,1,by=0.15); names(seq.ranges) <- as.character(1:length(seq.ranges))

foo <- to.plot %>%
  # .[score<=0] %>% 
  .[!is.na(motif_score)] %>%
  # merge(chip_background.dt,by="chip_anno") %>%
  .[,log_value:=log(value+1)] %>%
  # .[,log_value:=log(value+1)+abs(rnorm(.N,mean=0,sd=0.03))] %>% # add some noise because otherwise they are all zeros
  .[,c("chip_anno","idx","log_value","correlation_score","max_accessibility_score","motif_score","score")]


bar1 <-  foo %>%
  .[score<0] %>% .[,score:=abs(score)] %>%
  .[,score:=minmax.normalisation(score), by="chip_anno"] %>%
  .[,score_group:=cut(abs(score), breaks=seq.ranges)] %>% .[,score_group:=seq.ranges[as.numeric(score_group)]] %>%
  # Prepare for plotting
  melt(id.vars=c("chip_anno","idx","log_value"), measure.vars=c("score_group"), variable.name="model", value.name="predicted_value") %>%
  .[,N:=.N,by=c("model","chip_anno","predicted_value")] %>%
  .[,.(log_value=mean(log_value,na.rm=T), N=.N), by=c("model","chip_anno","predicted_value")] %>%
  # Filter settings with a small number of observations
  .[!is.na(predicted_value) & N>=3] %>% .[,class:="Negative"]

bar2 <-  foo %>%
  .[score>0] %>% .[,score:=abs(score)] %>%
  .[,score:=minmax.normalisation(score), by="chip_anno"] %>%
  .[,score_group:=cut(abs(score), breaks=seq.ranges)] %>% .[,score_group:=seq.ranges[as.numeric(score_group)]] %>%
  # Prepare for plotting
  melt(id.vars=c("chip_anno","idx","log_value"), measure.vars=c("score_group"), variable.name="model", value.name="predicted_value") %>%
  .[,N:=.N,by=c("model","chip_anno","predicted_value")] %>%
  .[,.(log_value=mean(log_value,na.rm=T), N=.N), by=c("model","chip_anno","predicted_value")] %>%
  # Filter settings with a small number of observations
  .[!is.na(predicted_value) & N>=3] %>% .[,class:="Positive"]

baz <-  rbind(bar1,bar2)

p <- ggline(baz, x="predicted_value", y="log_value", color="class", plot_type="b") +
  facet_wrap(~chip_anno, scales="free") +
  scale_color_brewer(palette="Set1") +
  # coord_cartesian(ylim=c(0,1)) +
  labs(y="Observed ChIP-seq signal", x="In silico binding score") +
  theme(
    strip.background = element_rect(colour="black", fill=NA),
    axis.text = element_text(size=rel(0.55)),
    axis.title = element_text(size=rel(0.8)),
    legend.title = element_blank(),
    legend.position = "top"
  )

pdf(file.path(io$outdir,"insilico_vs_chipseq_lineplot_negative_values.pdf"), width=7, height=5)
print(p)
dev.off()

##########
## Test ##
##########

# to.plot.compare_models[chip_anno=="TBX5"]
# unique(foo[chip_anno=="TBX5",motif_score])

# corr.peaks <- which(dropNA2matrix(assay(tf2peak_cor.se[,"POU5F1"],"cor"))[,1]>=opts$min.corr) %>% names

# library(BSgenome.Mmusculus.UCSC.mm10)
# BSgenome()
# BSgenome <- validBSgenome("BSgenome.Mmusculus.UCSC.mm10")
# seqlevels(subject) <- paste0("chr",levels(seqnames(subject)))
# 
# subject <- data.frame(chr="chr1", start=68590325, end=68590925) %>% makeGRangesFromDataFrame
# pwms <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/motif_seqlogo/PWMatrixList.rds")
# motifmatchr::matchMotifs(
#   pwms = pwms["TAL1"],
#   subject = subject,
#   genome = BSgenome.Mmusculus.UCSC.mm10, 
#   out = "positions", 
#   p.cutoff = 0.001, 
#   w = 7
# )
# 
# getSeq(Mmusculus, "chr1",68590647,68590656)
