here::i_am("atac/archR/differential/metacells/celltype/parse_differential_results.R")

# Load default settings
source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--diff_results_dir',   type="character",     help='File')
p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per group')
p$add_argument('--outdir',             type="character",     help='File')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$diff_results_dir <- file.path(io$basedir,"results/atac/archR/differential/metacells/celltype/PeakMatrix")
# args$min_cells <- 5
# args$outdir <- file.path(io$basedir,"results/atac/archR/differential/metacells/celltype/PeakMatrix/parsed")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

################################################
## Load differential expression and fetch TFs ##
################################################

stats.dt <- data.table(celltypeA=as.character(NA), celltypeB=as.character(NA), groupA_N=as.integer(NA), groupB_N=as.integer(NA), included=as.logical(NA))
diff_results_list <- list()

# i <- "Visceral_endoderm"; j <- "Surface_ectoderm"
for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    
    if (i!=j) {
      file <- file.path(args$diff_results_dir,sprintf("%s_vs_%s.txt.gz",opts$celltypes[[i]],opts$celltypes[[j]]))
      if (file.exists(file)) {
        tmp <- fread(file) %>% 
          .[,c("celltypeA","celltypeB"):=list(opts$celltypes[[i]],opts$celltypes[[j]])]
        
        if (tmp[celltypeA==opts$celltypes[[i]],groupA_N][1]>=args$min_cells & tmp[celltypeB==opts$celltypes[[j]],groupB_N][1]>=args$min_cells) {
          stats.dt <- rbind(stats.dt,data.table(celltypeA=opts$celltypes[[i]], celltypeB=opts$celltypes[[j]], groupA_N=tmp[celltypeA==opts$celltypes[[i]],groupA_N][1], groupB_N=tmp[celltypeB==opts$celltypes[[j]],groupB_N][1], included=TRUE))
          diff_results_list[[sprintf("%s_vs_%s",opts$celltypes[[i]],opts$celltypes[[j]])]] <- tmp
        } else {
          stats.dt <- rbind(stats.dt,data.table(celltypeA=opts$celltypes[[i]], celltypeB=opts$celltypes[[j]], groupA_N=tmp[celltypeA==opts$celltypes[[i]],groupA_N][1], groupB_N=tmp[celltypeB==opts$celltypes[[j]],groupB_N][1], included=FALSE)) 
        }
      } else {
        print(sprintf("%s not found...",file))
      }
    }
  }
}
 
##########
## Save ##
##########

fwrite(stats.dt[-1], file.path(args$outdir,"diff_stats.txt.gz"), sep="\t", quote=F, na="NA")
fwrite(rbindlist(diff_results_list), file.path(args$outdir,"diff_results.txt.gz"), sep="\t", quote=F, na="NA")

