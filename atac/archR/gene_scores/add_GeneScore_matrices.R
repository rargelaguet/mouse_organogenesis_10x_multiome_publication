here::i_am("atac/archR/gene_scores/add_GeneScore_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
# p$add_argument('--outdir',    type="character",    help='Output directory')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$metadata <- file.path(io$archR.directory,"sample_metadata_after_archR.txt.gz")
# args$outdir <- file.path(io$basedir,"results/atac/archR/gene_scores")
# args$threads <- 1
## END TEST ##

# dir.create(args$outdir, showWarnings = F)

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)

#####################
## Define gene set ##
#####################

genes.gr <- getGenes(ArchRProject)

#########################################
## Add Gene Scores using default model ##
#########################################

# Note that this will add the matrices to the arrowFiles
addGeneScoreMatrix(
  input = ArchRProject,
  genes = genes.gr,
  geneModel = "exp(-abs(x)/5000) + exp(-1)", # string should be a function of x, where x is the distance from the TSS.
  matrixName = "GeneScoreMatrix_distal",
  extendUpstream = c(1000, 1e+05),
  extendDownstream = c(1000, 1e+05),  # The minimum and maximum number of bp downstream of the transcription termination site to consider for gene activity score calculation.
  geneUpstream = 5000,                # Number of bp upstream the gene to extend the gene body.
  geneDownstream = 0,                 # Number of bp downstream the gene to extend the gene body.
  tileSize = 500,                     # The size of the tiles used for binning counts prior to gene activity score calculation.
  geneScaleFactor = 5,                # A numeric scaling factor to weight genes based on the inverse of their length 
  scaleTo = 10000,                    # Each column in the calculated gene score matrix will be normalized
  excludeChr = c("chrY", "chrM"),
  force = TRUE
)

# Save
# GeneScoreMatrix.se <- getMatrixFromProject(ArchRProject, binarize = FALSE, useMatrix = "GeneScoreMatrix_distal")
# rownames(GeneScoreMatrix.se) <- rowData(GeneScoreMatrix.se)$name
# saveRDS(GeneScoreMatrix.se, file.path(args$outdir,"GeneScoreMatrix_distal.rds"))

#############################################
## Add Gene Scores ignoring distal regions ##
#############################################


# TSS
addGeneScoreMatrix(
  input = ArchRProject,
  genes = genes.gr,
  useTSS = TRUE,
  extendTSS = TRUE,
  geneModel = "1",                    # string should be a function of x, where x is the distance from the TSS.
  matrixName = "GeneScoreMatrix_TSS",
  extendUpstream = c(0,0),
  extendDownstream = c(0,0),          # The minimum and maximum number of bp downstream of the transcription termination site to consider for gene activity score calculation.
  geneUpstream = 500,                # Number of bp upstream the gene to extend the gene body.
  geneDownstream = 100,                 # Number of bp downstream the gene to extend the gene body.
  tileSize = 100,                     # The size of the tiles used for binning counts prior to gene activity score calculation.
  geneScaleFactor = 1,                # A numeric scaling factor to weight genes based on the inverse of their length 
  scaleTo = 10000,                    # Each column in the calculated gene score matrix will be normalized
  excludeChr = c("chrY", "chrM"),
  force = TRUE
)

# Save
# GeneScoreMatrix.se <- getMatrixFromProject(ArchRProject, binarize = FALSE, useMatrix = "GeneScoreMatrix_TSS")
# rownames(GeneScoreMatrix.se) <- rowData(GeneScoreMatrix.se)$name
# saveRDS(GeneScoreMatrix.se, file.path(args$outdir,"GeneScoreMatrix_tss.rds"))

# Create a completion token
file.create(file.path(io["basedir"],"processed/atac/archR/addGeneScoreMatrix_completed.txt"))