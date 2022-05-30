library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$outfile <- file.path(io$basedir,"results_new/atac/archR/peak_calling/cpg_density_peaks.txt.gz")

################
## Load peaks ##
################

peaks.dt <- fread(io$archR.peak.metadata) %>%
  .[,c("chr","start","end")] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]# %>%
  # .[,c("idx")]


####################################
## Calculate CpG density per peak ##
####################################

chr_lengths.dt <- data.table(
  chr = unique(peaks.dt$chr),
  chr_length = seqlengths(Mmusculus) %>% .[unique(peaks.dt$chr)]
)
peaks.dt <- merge(peaks.dt, chr_lengths.dt, by="chr")

# Filter features that exceed chr  length
peaks.dt <- peaks.dt[end<chr_length]

# Get sequence
seq <- getSeq(Mmusculus, peaks.dt$chr, peaks.dt$start, peaks.dt$end+1)

# Calculate CpG density
peaks.dt$cpg_density <- round(dinucleotideFrequency(seq)[,"CG"]/width(seq),4)

##################
## Save results ##
##################

fwrite(peaks.dt[,c("idx","cpg_density")], io$outfile, col.names=T, quote=F, sep="\t")
