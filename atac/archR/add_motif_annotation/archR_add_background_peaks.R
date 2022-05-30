
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/add_motif_annotation/archR_add_background_peaks.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--method',        type="character", default="chromVAR",                              help='ArchR or chromVAR')
p$add_argument('--number_background_peaks',     type="integer",    default=50,    help='Number of background peaks')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$number_background_peaks <- 50
args$method <- "chromVAR"
args$threads <- 1
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))

########################
## Load ArchR Project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads = args$threads)

##########################
## Add background peaks ##
##########################

# This function will compute background peaks controlling for total accessibility and GC-content
# changes in the ArchR project: (1) it creates Background-Peaks.rds and (2) adds "bgdPeaks" entry to "metadata(getPeakSet(ArchRProject))"

# Background peaks are chosen by sampling peaks based on similarity in GC content and # of fragments across samples using the Mahalanobis distance. 
# - The "w" paramter controls how similar background peaks should be. 
# - The "binSize" parameter controls the precision with which the similarity is computed. Increasing "binSize" will make the function run slower.
# Returns a matrix with one row per peak and one column per iteration. values in a row represent indices of background peaks for the peak with that index

ArchRProject <- addBgdPeaks(
  ArchRProj = ArchRProject,
  nIterations = args$number_background_peaks,
  w = 0.1,
  binSize = 50,
  method = args$method,
  seed = 42,
  outFile = file.path(getOutputDirectory(ArchRProject), "Background-Peaks.rds"),   # default
  force = TRUE
)

# if (!file.exists(metadata(ArchRProject@peakSet)$bgdPeaks)) {

##########
## TEST ##
##########

R.utils::sourceDirectory("/bi/group/reik/ricard/scripts/git/archR/R", verbose=T, modifiedOnly=FALSE)

ArchRProj = ArchRProject
nIterations = 50
w = 0.1
binSize = 50
seed = 1
method = "chromVAR"

