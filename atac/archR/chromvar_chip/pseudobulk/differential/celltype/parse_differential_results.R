here::i_am("atac/archR/differential/pseudobulk/celltype/parse_differential_results.R")

# Load default settings
source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--diff_results_dir',   type="character",     help='File')
p$add_argument('--outfile',             type="character",     help='File')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$diff_results_dir <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltypes/CISBP")
# args$outfile <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltypes/CISBP/diff_results.txt.gz")
## END TEST ##

# I/O
dir.create(dirname(args$outfile), showWarnings=F, recursive=T)

################################################
## Load differential expression and fetch TFs ##
################################################

diff_results_list <- list()

# i <- "Visceral_endoderm"; j <- "Surface_ectoderm"
for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    
    if (i!=j) {
      file <- file.path(args$diff_results_dir,sprintf("%s_vs_%s.txt.gz",opts$celltypes[[i]],opts$celltypes[[j]]))
      if (file.exists(file)) {
        tmp <- fread(file) %>% .[,c("celltypeA","celltypeB"):=list(opts$celltypes[[i]],opts$celltypes[[j]])]
        diff_results_list[[sprintf("%s_vs_%s",opts$celltypes[[i]],opts$celltypes[[j]])]] <- tmp
      } else {
        print(sprintf("%s not found...",file))
      }
    }
  }
}
 
##########
## Save ##
##########

fwrite(rbindlist(diff_results_list), args$outfile, sep="\t", quote=F, na="NA")
