here::i_am("atac/archR/differential/pseudobulk/celltype_genotype/parse_differential_results.R")

# Load default settings
source(here::here("settings.R"))
# source(here::here("utils.R"))

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
# args$diff_results_dir <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype_genotype/PeakMatrix")
# args$outfile <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype_genotype/PeakMatrix/parsed/diff_results.txt.gz")
## END TEST ##

# I/O
dir.create(dirname(args$outfile), showWarnings = F, recursive = T)

##########################################
## Load differential expression results ##
##########################################

diff_results_list <- list()

# i <- "Visceral_endoderm"; j <- "Surface_ectoderm"
for (i in 1:length(opts$celltypes)) {
  file <- file.path(args$diff_results_dir,sprintf("%s.txt.gz",opts$celltypes[[i]]))
  if (file.exists(file)) {
    tmp <- fread(file) %>% .[,celltype:=opts$celltypes[[i]]]
    if (nrow(tmp)>1) {
      diff_results_list[[opts$celltypes[[i]]]] <- tmp
    }
  } else {
    print(sprintf("%s not found...",file))
  }
}
 
print(names(diff_results_list))

##########
## Save ##
##########

fwrite(rbindlist(diff_results_list), args$outfile, sep="\t", quote=F, na="NA")

