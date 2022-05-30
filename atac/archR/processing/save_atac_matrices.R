here::i_am("atac/archR/processing/save_atac_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--matrix',      type="character",  help='Matrix to save')
p$add_argument('--outfile',     type="character",  help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$matrix <- "GeneScoreMatrix_TSS"
# args$outfile <- file.path(io$basedir,sprintf("processed/atac/archR/Matrices/%s_summarized_experiment.rds",args$matrix))
## END TEST ##

print(args)

# I/O
dir.create(dirname(args$outfile), showWarnings=F)

###################
## Load metadata ##
###################

cells_metadata.dt <- fread(args$metadata) %>%
  # .[pass_atacQC==TRUE & doublet_call==FALSE]
  .[pass_atacQC==TRUE]

########################
## Load ArchR Project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = 1)

ArchRProject <- loadArchRProject(args$archr_directory)[cells_metadata.dt$cell]

# Sanity checks
# mean(rownames(ArchRProject)%in%cells_metadata.dt$cell)
# mean(cells_metadata.dt$cell%in%rownames(ArchRProject))
# table(cells_metadata.dt[!cell%in%rownames(ArchRProject),sample])
stopifnot(args$matrix %in% getAvailableMatrices(ArchRProject))

################
## PeakMatrix ##
################

if (args$matrix=="PeakMatrix") {

	atac.se <- getMatrixFromProject(ArchRProject, binarize = FALSE, useMatrix = "PeakMatrix")

	# Define peak names
	row_ranges.dt <- rowRanges(atac.se) %>% as.data.table %>% 
	  setnames("seqnames","chr") %>%
	  .[,c("chr","start","end")] %>%
	  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
	rownames(atac.se) <- row_ranges.dt$idx

}

#####################
## GeneScoreMatrix ##
#####################

if (grepl("GeneScoreMatrix",args$matrix)) {

	atac.se <- getMatrixFromProject(ArchRProject, binarize = FALSE, useMatrix = args$matrix)

	# Define gene names
	rownames(atac.se) <- rowData(atac.se)$name

	# Filter genes
	# atac.se <- atac.se[grep("^Rik|Rik$|^mt-|^Rps-|^Rpl-|^Gm|^Mir|^Olfr",rownames(atac.se),invert=T),]
}

##########
## Save ##
##########

# Sanity checks
stopifnot(sum(duplicated(rownames(atac.se)))==0)

saveRDS(atac.se, args$outfile)