source(here::here("settings.R"))

library(csaw)
library(GenomicRanges)
library(BiocParallel)

#####################
## Define settings ##
#####################

# Options
# opts$samples <- c(
#   "p8_Control_REP1",
#   "p8_Control_REP2",
#   "p8_Control_REP3",
#   "p8_TG_REP1",
#   "p8_TG_REP2",
#   "p8_TG_REP3",
#   "p11_Control_REP1",
#   "p11_Control_REP2",
#   "p11_Control_REP3",
#   "p11_TG_REP1",
#   "p11_TG_REP2",
#   "p11_TG_REP3",
#   "p14_Control_REP1",
#   "p14_Control_REP2",
#   "p14_Control_REP3",
#   "p14_TG_REP1",
#   "p14_TG_REP2",
#   "p14_TG_REP3",
#   "p19_Control_REP1",
#   "p19_Control_REP2",
#   "p19_Control_REP3",
#   "p19_TG_REP1",
#   "p19_TG_REP2",
#   "p19_TG_REP3",
#   "p22_Control_REP1",
#   "p22_Control_REP2",
#   "p22_Control_REP3",
#   "p22_TG_REP1",
#   "p22_TG_REP2",
#   "p22_TG_REP3",
#   "p24_Control_REP1",
#   "p24_Control_REP2",
#   "p24_Control_REP3",
#   "p24_TG_REP1",
#   "p24_TG_REP2",
#   "p24_TG_REP3",
#   "p29_Control_REP1",
#   "p29_Control_REP2",
#   "p29_Control_REP3",
#   "p29_TG_REP1",
#   "p29_TG_REP2",
#   "p29_TG_REP3",
#   "pA_Control_REP1",
#   "pA_Control_REP2",
#   "pA_Control_REP3",
#   "pA_TG_REP1",
#   "pA_TG_REP2",
#   "pA_TG_REP3"
# )

opts$samples <- c(
  "p8_Control", 
  "p11_Control", 
  "p14_Control", 
  "p19_Control", 
  "p22_Control", 
  "p24_Control", 
  "p29_Control", 
  "p8_TG", 
  "p11_TG", 
  "p14_TG", 
  "p19_TG", 
  "p22_TG", 
  "p24_TG", 
  "p29_TG"
)

# opts$chr <- paste0("chr",c(1:22,"X","Y"))
opts$chr <- c(1:22,"X","Y")
opts$window_width <- 5000
opts$window_spacing <- 2500

# Parallel processing
BPPARAM <- BiocParallel::MulticoreParam(workers=round(length(opts$samples)/6))

# I/O
io$bam_files <- file.path(io$basedir,sprintf("bam/atac/merged/%s.bam",opts$samples))
io$blacklist <- "/share/rargelaguet/data/hg38_regulation/blacklist/hg38-blacklist.v2.bed"
stopifnot(file.exists(io$bam_files))

io$outdir <- file.path(io$basedir,"processed/atac/csaw"); dir.create(io$outdir, showWarnings=F, recursive=T)


############################
## Load Blacklist regions ##
############################

blacklist <- read.csv(io$blacklist, sep = "\t", header = F)
blacklist.gr <- GenomicRanges::GRanges(seqnames = blacklist$V1, ranges = IRanges(start = blacklist$V2 + 1, end = blacklist$V3))

##########
## csaw ##
##########

# list bam files
# work_dir <- "/home/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/bam"
# list.bam <- list.files(path = work_dir, pattern = ".bam$", full.names = T, recursive = T)
# list.bam <- list.bam[gsub(pattern = ".bam", replacement = "", x = basename(list.files(path = work_dir, pattern = ".bam$", full.names = T, recursive = T))) %in% samples]

# Examine the distribution of fragment sizes
# out <- csaw::getPESizes(bam.file = list.bam[1])
# frag.sizes <- out$sizes[out$sizes <= 800]
# hist(frag.sizes, breaks = 50, xlab = "Fragment sizes (bp)", ylab = "Frequency", main = "", col = "grey80")

# Construct a readParam object to standardize the parameter settings in this analysis
# Set the minimum mapping quality score to 10 to remove poorly or non-uniquely aligned reads. 
# Paired-end datasets are accomodated by setting pe=“both” in the param object.
param <- csaw::readParam(minq = 10, discard = blacklist.gr, restrict = opts$chr, pe = "both", max.frag = 500)

# Count backgroud bins to calculate normalization factors:
# bins.10kb <- csaw::windowCounts(list.bam, bin = TRUE, width = 10000, param = param)
# bins.10kb <- csaw::normFactors(bins.10kb, se.out = bins.10kb)

# Count reads over sliding windows
# atac_windows.rse <- csaw::windowCounts(bam.files = io$bam_files, param = param, width = opts$window_width, spacing = opts$window_spacing, ext = 100, BPPARAM = BPPARAM)
atac_windows.rse <- csaw::windowCounts(bam.files = io$bam_files, param = param, width = opts$window_width, spacing = opts$window_spacing, ext = 100)
atac_windows.rse <- csaw::normFactors(atac_windows.rse)

# TO DO: ADD ROWNAMES AND COLNAMES

# Save
saveRDS(atac_windows.rse, file.path(io$outdir,sprintf("windows_width%s_spacing%s.rds",opts$window_width,opts$window_spacing)))
