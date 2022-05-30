here::i_am("atac/archR/add_motif_annotation/plot_motif_seqlogo.R")

source(here::here("settings.R"))

suppressPackageStartupMessages(library(TFBSTools))
suppressPackageStartupMessages(library(ggseqlogo))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',             type="character",            help='Motif annotation')
# p$add_argument('--peak_annotation_file',          type="character",                               help='')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--test',    action="store_true",             help='Test mode')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$motif_annotation <- "CISBP"
args$peak_annotation_file <- file.path(io$archR.directory,"Annotations/peakAnnotation.rds")
args$outdir <- file.path(io$archR.directory,sprintf("Annotations/seqlogo/%s",args$motif_annotation))
args$test <- FALSE
## END TEST ##

# Parse arguments
dir.create(args$outdir, showWarnings=F, recursive = T)

###############
## Load data ##
###############

# raw position frequency matrix (PFM)
pwm <- readRDS(args$peak_annotation_file)[[args$motif_annotation]][["motifs"]]
# pwm <- readRDS(io$peak_annotation_file)[["Motif_JASPAR2020_human"]][["motifs"]]

##########
## Plot ##
##########

motifs.to.plot <- names(pwm)
if (args$test) {
  motifs.to.plot <- motifs.to.plot %>% head(n=3)
}

# postProbs = (PFM + bg * pseudocounts) / (colSums(PFM) + sum(bg) * pseudocounts)
# priorProbs = bg / sum(bg)
# PWM_log2probratio = log2(postProbs / priorProbs)

grep("YBX",motifs.to.plot,value=T)

# i <- "YBX2_827"
for (i in motifs.to.plot) {

  # position weight matrix (PWM)
  
  if (args$motif_annotation=="JASPAR") {
    tmp <- toPWM(pwm[[i]], type="prob") %>% as.matrix
  } else if (args$motif_annotation=="CISBP") {
    tmp <- (2**as.matrix(pwm[[i]]))*0.25  # this is not entirely accurate
  }
  
  p <- ggseqlogo(tmp) + 
    theme(
      axis.line = element_line(size=rel(0.5), color="black"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=rel(0.75)),
      axis.title.y = element_text(size=rel(0.75)),
      # axis.title.y = element_blank()
    )
  pdf(file.path(args$outdir,sprintf("seqlogo_%s_%s.pdf",args$motif_annotation,i)), width=5, height=2.2)
  print(p)
  dev.off()
}

# Completion token
file.create(file.path(args$outdir,"completed.txt"))