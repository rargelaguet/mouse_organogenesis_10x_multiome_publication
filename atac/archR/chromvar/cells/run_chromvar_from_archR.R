# https://www.archrproject.com/bookdown/chromvar-deviatons-enrichment-with-archr.html
# Note: this requires the previous execution of (...)/add_motif_annotation/archR_add_motif_annotation.R
here::i_am("atac/archR/chromvar/cells/run_chromvar_from_archR.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(GenomicRanges))
suppressMessages(library(ArchR))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--motif_annotation',  type="character",                help='Motif annotation') 
p$add_argument('--metadata',  type="character",                help='Metadata') 
p$add_argument('--threads',            type="integer",    default=1,    help='Number of cores')
p$add_argument('--force', action="store_true", 				help='Force')
p$add_argument('--outdir',  type="character",                help='Output directory') 
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$motif_annotation <- "CISBP"
# args$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$threads <- 1
# args$force <- FALSE
# args$outdir <- file.path(io$basedir,"results/atac/archR/chromvar")
## END TEST ##

#####################
## Load metadata ##
#####################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]

########################
## Load ArchR Project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]

#######################################
## Sanity checks on the ArchR object ##
#######################################

# Load Peak Set
# if (file.exists(io$archR.peakSet.granges)) {
# 	ArchRProject <- addPeakSet(ArchRProject, peakSet = readRDS(io$archR.peakSet.granges), force = TRUE)
# }

# Load motif annotations over peaks
tmp <- file.path(args$archr_directory,"Annotations/peakAnnotation.rds")
if (file.exists(tmp)) {
	ArchRProject@peakAnnotation <- readRDS(tmp)
}
stopifnot(args$motif_annotation%in%names(ArchRProject@peakAnnotation))

##########################
## Add background peaks ##
##########################

# This function will compute background peaks controlling for total accessibility and GC-content
# changes in the ArchR project: (1) it creates Background-Peaks.rds and (2) adds "bgdPeaks" entry to "metadata(getPeakSet(ArchRProject))"

# Background peaks are chosen by sampling peaks based on similarity in GC content and # of fragments across samples using the Mahalanobis distance. 
# The "w" paramter controls how similar background peaks should be. The bs parameter controls the precision with which the similarity is computed; 
# increasing bs will make the function run slower.
# Returns a matrix with one row per peak and one column per iteration. values in a row represent indices of background peaks for the peak with that index

# if (!file.exists(metadata(ArchRProject@peakSet)$bgdPeaks)) {
#   ArchRProject <- addBgdPeaks(
#     ArchRProject,
#     nIterations = 50,
#     w = 0.1,
#     binSize = 50,
#     method = "chromVAR",
#     seed = 42,
#     # outFile = file.path(getOutputDirectory(ArchRProj), "Background-Peaks.rds"),
#     force = TRUE
#   )
# }

metadata(ArchRProject@peakSet)$bgdPeaks <- file.path(args$archr_directory, "Background-Peaks.rds")

stopifnot("bgdPeaks"%in%names(metadata(ArchRProject@peakSet)))

###################################
## Compute deviations for motifs ##
###################################

# The function computeDeviations returns a SummarizedExperiment with two "assays":  
# - The first matrix (accessible via `deviations(dev)` or `assays(dev)$deviations)` will give the bias corrected deviation in accessibility for each set of peaks (rows) for each cell or sample (columns). This metric represent how accessible the set of peaks is relative to the expectation based on equal chromatin accessibility profiles across cells/samples, normalized by a set of background peak sets matched for GC and average accessibility. 
# - The second matrix `deviationScores(dev)` or `assays(deviations)$z` gives the deviation Z-score, which takes into account how likely such a score would occur if randomly sampling sets of beaks with similar GC content and average accessibility.

print("Running chromVAR implementation in ArchR...")

matrix_name <- paste0("DeviationMatrix_",args$motif_annotation)

if (matrix_name %in% getAvailableMatrices(ArchRProject)) {
	if (args$force) {
		ArchRProject <- addDeviationsMatrix(
		  ArchRProject, 
		  matrixName = matrix_name,
		  peakAnnotation = args$motif_annotation,
		  out = c("z", "deviations"),
		  binarize = FALSE,
		  force = TRUE
		)
	} else {
		stop(sprintf("%s already found in the ArchR object. Pass the --force argument to replace the existing data matrix...",matrix_name))
	}

} else {
	ArchRProject <- addDeviationsMatrix(
	  ArchRProject, 
	  matrixName = matrix_name,
	  peakAnnotation = args$motif_annotation,
	  out = c("z", "deviations"),
	  binarize = FALSE
	)
}

##########
## Save ##
##########

# Fetch deviations matrix as a SummarisedExperiment
chromvar_deviations.se <- getMatrixFromProject(ArchRProject, matrix_name)

# Save SummarisedExperiment object
saveRDS(chromvar_deviations.se, file.path(args$outdir,sprintf("chromVAR_deviations_%s_archr.rds",args$motif_annotation)))

##########
## TEST ##
##########

# if (grepl("ricard",Sys.info()['nodename'])) {
#   R.utils::sourceDirectory("/Users/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)
# } else if (grepl("ebi",Sys.info()['nodename'])) {
#   R.utils::sourceDirectory("/homes/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)
# } else {
#   stop("Computer not recognised")
# }

# ArchRProj = ArchRProject
# peakAnnotation = "Motif_JASPAR2020_human"
# matches = NULL
# bgdPeaks = getBgdPeaks(ArchRProj, method = "chromVAR")
# matrixName = NULL
# out = c("z", "deviations")
# binarize = FALSE
# threads = 1
# verbose = TRUE
# parallelParam = NULL
# force = FALSE
# logFile = createLogFile("addDeviationsMatrix")
