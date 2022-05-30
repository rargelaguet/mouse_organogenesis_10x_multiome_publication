here::i_am("rna/differential/metacells/celltype/run_diff_expr_celltype.R")

# Load default settings
source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",     help='SingleCellExperiment file')
p$add_argument('--outdir',          type="character",     help='Output directory')
p$add_argument('--test_mode',       action="store_true",  help='Test mode? subset data')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk_with_replicates.rds")
# args$outdir <- file.path(io$basedir,"results/rna/differential/pseudobulk/with_replicates/celltype")
# args$test_mode <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

io$script <- here::here("rna/differential/pseudobulk/with_replicates/celltype/differential_rna_celltype_pseudobulk.R")
dir.create(args$outdir, showWarnings=FALSE, recursive=TRUE)

#########
## Run ##
#########

if (args$test_mode) {
  print("Test mode activated, running only a few comparisons...")
  opts$celltypes <- opts$celltypes %>% head(n=3)
}

for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      groupA <- opts$celltypes[[i]]
      groupB <- opts$celltypes[[j]]
      
      outfile <- sprintf("%s/%s_vs_%s.txt.gz", args$outdir,groupA,groupB)

      # Define LSF command
      if (grepl("BI",Sys.info()['nodename'])) {
        lsf <- ""
      } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
        lsf <- sprintf("sbatch -n 1 --mem 7G --wrap")
      }
      cmd <- sprintf("%s 'Rscript %s --sce %s --groupA %s --groupB %s --outfile %s'", 
        lsf, io$script, args$sce, groupA, groupB, outfile)
      # if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")

      # Run
      print(cmd)
      system(cmd)
    }
  }
}


# Completion token
file.create(file.path(args$outdir,"completed.txt"))