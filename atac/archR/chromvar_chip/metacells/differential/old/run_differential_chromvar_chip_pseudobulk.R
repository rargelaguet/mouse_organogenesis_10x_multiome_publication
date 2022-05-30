suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/chromvar_chip/pseudobulk/differential/run_differential_chromvar_chip_pseudobulk.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--chromvar_chip_pseudobulk',  type="character",              help='Motif annotation') 
p$add_argument('--outdir',  type="character",              help='Motif annotation') 
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args$motif_annotation <- "JASPAR"
# args$chromvar_chip_pseudobulk <- file.path(io$basedir,sprintf("results_new/atac/archR/chromvar_chip/pseudobulk/chromVAR_deviations_%s_archr_chip.rds",args$motif_annotation))
# args$outdir <- file.path(io$basedir,sprintf("results_new/atac/archR/chromvar_chip/pseudobulk/differential/%s",args$motif_annotation))
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)

#####################################
## Load pseudobulk chromVAR scores ##
#####################################

chromvar_deviations_pseudobulk.se <- readRDS(args$chromvar_chip_pseudobulk)
  
######################################
## Differential motif accessibility ##
######################################

# i <- 1; j <- 2
for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      foo <- assay(chromvar_deviations_pseudobulk.se[,opts$celltypes[[j]]])[,1]
      bar <- assay(chromvar_deviations_pseudobulk.se[,opts$celltypes[[i]]])[,1]
      
      chromvar_diff.dt <- data.table(
        gene = names(foo), 
        diff = round(foo-bar,2) 
        # groupA = opts$celltypes[[i]], 
        # groupB = opts$celltypes[[j]]
      ) %>% sort.abs("diff") 
      
      # save      
      outfile <- file.path(args$outdir,sprintf("%s_vs_%s_%s_chromVAR_chip_pseudobulk.txt.gz",opts$celltypes[[i]],opts$celltypes[[j]],args$motif_annotation))
      fwrite(chromvar_diff.dt, outfile, sep="\t")
    }
  }
}
