here::i_am("atac/archR/differential/metacells/celltype/run_diff_acc_celltype_metacells.R")

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--metadata',          type="character",   help='')
p$add_argument('--atac_matrix_file',        type="character",                               help='')
p$add_argument('--group_variable',          type="character",   help='')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--min_cells',       type="integer",       default=50,      help='Minimum number of cells per cell type')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset data')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$matrix <- "PeakMatrix" # "GeneScoreMatrix_TSS"
# args$metadata <- file.path(io$basedir,"results/atac/archR/metacells/metacells_metadata.txt.gz")
# args$atac_matrix_file <- file.path(io$basedir,sprintf("results/atac/archR/metacells/%s_summarized_experiment_metacells.rds",args$matrix))
# args$group_variable <- "celltype"
# args$min_cells <- 5
# args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/metacells/%s/%s",args$group_variable,args$matrix))
# args$test_mode <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
io$script <- here::here("atac/archR/differential/metacells/run_differential_accessibility_metacells.R")
dir.create(args$outdir, showWarnings=FALSE, recursive=TRUE)

########################
## Load cell metadata ##
########################

metacells_metadata.dt <- fread(args$metadata) %>%
  .[,celltype.mapped:=celltype]

stopifnot(args$group_variable%in%colnames(metacells_metadata.dt))

# subset celltypes with sufficient number of cells
metacells_metadata.dt <- fread(args$metadata) %>% 
  setnames(args$group_variable,"celltype") %>%
  .[,N:=.N,by="celltype"] %>% .[N>=args$min_cells] %>% .[,N:=NULL]
celltypes.to.use <- unique(metacells_metadata.dt$celltype)# %>% head(n=3)

# print stats
celltype.stats <- table(metacells_metadata.dt$celltype)
print(celltype.stats)

#########
## Run ##
#########

if (args$test_mode) {
  print("Test mode activated, running only a few comparisons...")
  celltypes.to.use <- celltypes.to.use %>% head(n=3)
}

stats.dt <- data.table(groupA=as.character(NA), groupB=as.character(NA), N_A=as.integer(NA), N_B=as.integer(NA))

for (i in 1:length(celltypes.to.use)) {
  for (j in i:length(celltypes.to.use)) {
    if (i!=j) {
      groupA <- celltypes.to.use[[i]]
      groupB <- celltypes.to.use[[j]]
      # print(sprintf("i=%s (%s), j=%s (%s)",i,groupA,j,groupB))
      
      outfile <- sprintf("%s/%s_%s_vs_%s.txt.gz", args$outdir,args$matrix,groupA,groupB)
      if (!file.exists(outfile)) {

        # Define LSF command
        if (grepl("BI",Sys.info()['nodename'])) {
          lsf <- ""
        } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
          lsf <- sprintf("sbatch -n 1 --mem 7G --wrap")
        }
        cmd <- sprintf("%s 'Rscript %s --atac_matrix_file %s --metadata %s --groupA %s --groupB %s --matrix %s --group_variable %s --outfile %s'", 
          lsf, io$script, args$atac_matrix_file, args$metadata, groupA, groupB, args$matrix, args$group_variable, outfile)

        # Run
        print(cmd)
        system(cmd)
        
        # save stats
        stats.dt <- rbind(stats.dt, data.table(groupA=groupA, groupB=groupB, N_A=celltype.stats[[groupA]], N_B=celltype.stats[[groupB]]))
        
      } else {
        print(sprintf("%s already exists...",outfile))
      }
    }
  }
}


# Save stats
fwrite(stats.dt[-1], file.path(args$outdir,"diff_stats.txt"), sep="\t", quote=F)

# Completion token
file.create(file.path(args$outdir,"completed.txt"))