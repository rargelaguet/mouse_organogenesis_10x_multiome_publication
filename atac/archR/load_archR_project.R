suppressPackageStartupMessages(library(ArchR))

#####################
## Define settings ##
#####################

io$archR.directory <- file.path(io$basedir,"processed/atac/archR")
# io$archR.directory <- file.path(io$basedir,"test/processed/atac/archR")
# io$atac.peak.annotation <- file.path(io$basedir,"/original/atac_peak_annotation.tsv")
io$archR.projectMetadata <- file.path(io$archR.directory,"projectMetadata.rds")
io$archR.peakSet.granges <- file.path(io$archR.directory,"PeakSet.rds")

setwd(io$archR.directory)

####################
## Define options ##
####################

addArchRGenome("mm10")
addArchRThreads(threads = 1) 

########################
## Load ArchR project ##
########################

ArchRProject <- loadArchRProject(io$archR.directory)

# Load ArchR projectMetadata
if (file.exists(io$archR.projectMetadata)) {
	ArchRProject@projectMetadata <- readRDS(io$archR.projectMetadata)
}

# Load peaks
if (file.exists(io$archR.peakSet.granges)) {
	ArchRProject <- addPeakSet(ArchRProject, peakSet = readRDS(io$archR.peakSet.granges), force = TRUE)
}

# Load motif annotations over peaks
if (!is.null(ArchRProject@peakAnnotation)) {
	io$archR.peakAnnotation <- file.path(io$archR.directory,"Annotations/peakAnnotation.rds")
	if (file.exists(io$archR.peakAnnotation)) {
		ArchRProject@peakAnnotation <- readRDS(io$archR.peakAnnotation)
	}
}

# Add background peaks
if (!is.null(getPeakSet(ArchRProject))) {
	io$archR.bgdPeaks <- file.path(io$archR.directory, "Background-Peaks.rds")
	if (!"bgdPeaks" %in% metadata(getPeakSet(ArchRProject))$bgdPeaks) {
		if (file.exists(io$archR.bgdPeaks)) metadata(ArchRProject@peakSet)$bgdPeaks <- io$archR.bgdPeaks
	}
}

##########
## TEST ##
##########

# ArchRProject@peakSet <- readRDS(io$archR.peakSet.granges)
# seqlevels(ArchRProject@peakSet) <- sort(seqlevels(ArchRProject@peakSet))
# ArchRProject@peakSet <- sort(ArchRProject@peakSet)


# getAvailableMatrices(ArchRProject)

# io$arrow.files <- opts$samples %>% 
#   # map_chr(~ sprintf("%s/%s.arrow",io$archR.directory,.))
#   map_chr(~ sprintf("%s.arrow",.))
# 
# ArchRProject <- ArchRProject(
#   ArrowFiles = io$arrow.files,
#   # outputDirectory = "ArchROutput",
#   outputDirectory = io$archR.directory,
#   copyArrows = FALSE
# )
# saveArchRProject(ArchRProject)
