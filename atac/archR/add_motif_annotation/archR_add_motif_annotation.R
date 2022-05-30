# https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/gene_scores/add_GeneScore_matrices.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
# p$add_argument('--outdir',    type="character",    help='Output directory')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
# args$threads <- 1
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

########################
## Load ArchR Project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads = args$threads)

##########################
## Add motif annotation ##
##########################

# cisbp (stringent threshold)
ArchRProject <- addMotifAnnotations(
  ArchRProject,
  motifSet = "cisbp",
  name = "Motif_cisbp",
  cutOff = 5e-05,
  width = 7,
  force = FALSE
)

# cisbp (lenient threshold)
# ArchRProject <- addMotifAnnotations(
#   ArchRProject,
#   motifSet = "cisbp",
#   name = "Motif_cisbp_lenient",
#   cutOff = 1e-04,
#   width = 7,
#   force = TRUE
# )

# homer
# ArchRProject <- addMotifAnnotations(
#   ArchRProject,
#   motifSet = "homer",
#   cutOff = opts$motif.pvalue.cutoff,
#   name = "Motif_homer",
#   force = TRUE
# )


# JASPAR2020 human (stringent)
ArchRProject <- addMotifAnnotations(
  ArchRProject, 
  motifSet = "JASPAR2020",      
  collection = "CORE",  
  species = "Homo sapiens",
  cutOff = 5e-05,   
  name = "Motif_JASPAR2020",
  force = FALSE
)

# JASPAR2020 human (lenient)
# ArchRProject <- addMotifAnnotations(
#   ArchRProject, 
#   motifSet = "JASPAR2020",      
#   collection = "CORE",  
#   species = "Homo sapiens",
#   cutOff = 1e-04,   
#   name = "Motif_JASPAR2020_lenient",
#   force = TRUE
# )


################################
## Save peakAnnotation object ##
################################

saveRDS(ArchRProject@peakAnnotation, sprintf("%s/Annotations/peakAnnotation.rds",io$archR.directory))
