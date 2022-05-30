here::i_am("atac/archR/differential/metacells/genotype/run_diff_acc_genotype_metacells.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',          type="character",   help='')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--atac_matrix_file',        type="character",                               help='')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--min_cells',       type="integer",       default=50,      help='Minimum number of cells per cell type')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset data')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$matrix <- "PeakMatrix" # "GeneScoreMatrix_TSS"
# args$metadata <- file.path(io$basedir,sprintf("results/atac/archR/metacells/all_cells/%s/metacells_metadata.txt.gz",args$matrix))
# args$atac_matrix_file <- file.path(io$basedir,sprintf("results/atac/archR/metacells/all_cells/%s/%s_summarized_experiment_metacells.rds",args$matrix,args$matrix))
# args$min_cells <- 5
# args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/metacells/genotype/%s",args$matrix))
# args$test_mode <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
io$script <- here::here("atac/archR/differential/metacells/run_differential_accessibility_metacells.R")
dir.create(args$outdir, showWarnings=FALSE, recursive=TRUE)

# Options
opts$samples <- c(
  "E8.5_CRISPR_T_KO",
  "E8.5_CRISPR_T_WT"
)

########################
## Load cell metadata ##
########################

metacells_metadata.dt <- fread(args$metadata) %>%
  .[sample%in%opts$samples & !is.na(genotype) & !is.na(celltype)] %>%
  .[,celltype_genotype:=sprintf("%s_%s",celltype,genotype)]

print(table(metacells_metadata.dt$sample))
print(table(metacells_metadata.dt$celltype))
print(table(metacells_metadata.dt$genotype))

######################
## Filter metacells ##
######################

# Only consider cell types with sufficient number of metacells
stats.dt <- metacells_metadata.dt[,.N,by=c("celltype","genotype")] %>% dcast(celltype~genotype, value.var="N", fill=0)
celltypes.to.use <- stats.dt[T_KO>=args$min_cells & WT>=args$min_cells,celltype] 
stats.dt <- stats.dt[celltype%in%celltypes.to.use]
print(stats.dt)

#########
## Run ##
#########

if (args$test_mode) {
  print("Test mode activated, running only a few comparisons...")
  celltypes.to.use <- celltypes.to.use %>% head(n=3)
}

# j <- "NMP"
for (j in celltypes.to.use) {
  outfile <- sprintf("%s/%s_WT_vs_KO.txt.gz", args$outdir,j); dir.create(dirname(outfile), showWarnings = F)
  if (!file.exists(outfile)) {
    
    # Define LSF command
    if (grepl("BI",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
      lsf <- sprintf("sbatch -n 1 --mem 8G --wrap")
    }
    cmd <- sprintf("%s 'Rscript %s --metadata %s --atac_matrix_file %s --samples %s --celltypes %s --groupA WT --groupB T_KO --matrix %s --group_variable genotype --outfile %s'",
      lsf, io$script, args$metadata, args$atac_matrix_file, paste(opts$samples,collapse=" "), j, args$matrix, outfile)
    
    # Run
    print(cmd)
    system(cmd)
  }
}


# Save stats
fwrite(stats.dt, file.path(args$outdir,"diff_stats.txt"), sep="\t", quote=F)

# Completion token
file.create(file.path(args$outdir,"completed.txt"))