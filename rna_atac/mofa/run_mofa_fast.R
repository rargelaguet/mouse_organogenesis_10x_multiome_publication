here::i_am("rna_atac/mofa/run_mofa_fast.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


suppressPackageStartupMessages(library(MOFA2))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--atac_dimred',          type="character",  help='')
p$add_argument('--rna_dimred',          type="character",   help='')
p$add_argument('--stages',       type="character",  default="all",  nargs='+',  help='Stages to plot')
p$add_argument('--samples',       type="character",  default="all",  nargs='+',  help='Samples to plot')
p$add_argument('--factors',           type="integer",    default=30,                  help='Number of MOFA factors')
p$add_argument('--remove_ExE_cells',       type="character",  default="False",  help='Remove ExE cells? ("True"/"False")')
p$add_argument('--outdir',          type="character",                               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# io$basedir <- file.path(io$basedir,"test")
# args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
# args$atac_dimred <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/PeakMatrix/remove_ExE_cells_False/batch_correction_None/lsi_nfeatures25000_ndims50.txt.gz")
# args$rna_dimred <- file.path(io$basedir,"results/rna/dimensionality_reduction/sce/batch_correction_by_sample_remove_ExE_cells_False/pca_features2500_pcs50.txt.gz")
# args$stages <- "all"
# args$samples <- "all"
# args$remove_ExE_cells <- "True"
# args$factors <- 30
# args$outdir <- file.path(io$basedir,"results/rna_atac/mofa/fast")
## END TEST ##

dir.create(args$outdir, showWarnings = F)

#####################
## Parse arguments ##
#####################

if (args$stages[1]=="all") {
  args$stages <- opts$stages
} else {
  stopifnot(args$stages%in%opts$stages)
}

if (args$remove_ExE_cells=="True") {
  args$remove_ExE_cells <- TRUE
} else if (args$remove_ExE_cells=="False") {
  args$remove_ExE_cells <- FALSE 
} else {
  stop('remove_ExE_cells should be "True" or "False"')
}

print(args)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  # .[(pass_atacQC==TRUE & pass_rnaQC==TRUE) & doublet_call==FALSE] %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[stage%in%args$stages]

opts$rna.cells <- sample_metadata[pass_rnaQC==TRUE,cell]
opts$atac.cells <- sample_metadata[pass_atacQC==TRUE,cell]

###############################################
## Load precomputed dimensionality reduction ##
###############################################

# RNA (PCA)
tmp <- fread(io$pca.rna)
opts$rna.cells <- intersect(tmp$cell,opts$rna.cells)
rna.mtx <- tmp %>% matrix.please %>% .[opts$rna.cells,] %>% t
rm(tmp)

# ATAC (LSI)
tmp <- fread(args$atac_dimred)
opts$atac.cells <- intersect(tmp$cell,opts$atac.cells)
atac.mtx <- tmp %>% matrix.please %>% .[opts$atac.cells,] %>% t
rm(tmp)

###########################
## Prepare data for MOFA ##
###########################

rna.mtx <- .augment_matrix(rna.mtx, unique(opts$rna.cells,opts$atac.cells))
atac.mtx <- .augment_matrix(atac.mtx, unique(opts$rna.cells,opts$atac.cells))

########################
## Create MOFA object ##
########################

MOFAobject <- create_mofa_from_matrix(list("RNA" = rna.mtx, "ATAC" = atac.mtx))
MOFAobject

####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$use_float32 <- TRUE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 30
model_opts$spikeslab_weights <- FALSE
model_opts$ard_weights <- FALSE

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
# train_opts$maxiter <- 5

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#####################
## Train the model ##
#####################

MOFAobject <- run_mofa(MOFAobject)

#########################
# Add samples metadata ##
#########################

metadata.to.mofa <- sample_metadata %>% copy %>%
  setnames("sample","batch") %>% setnames("cell","sample") %>%
  .[sample%in%unlist(samples_names(MOFAobject))] %>%
  setkey(sample) %>% .[unlist(samples_names(MOFAobject))]
samples_metadata(MOFAobject) <- metadata.to.mofa

##########
## Save ##
##########

saveRDS(MOFAobject, file.path(args$outdir,"mofa.rds"))
fwrite(metadata.to.mofa, file.path(args$outdir,"sample_metadata.txt.gz"), quote=F, sep="\t", na="NA")
