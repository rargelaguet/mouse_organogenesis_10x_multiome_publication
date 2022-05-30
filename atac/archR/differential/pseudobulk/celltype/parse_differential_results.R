here::i_am("atac/archR/differential/pseudobulk/celltype/parse_differential_results.R")

# Load default settings
source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--diff_results_dir',   type="character",     help='File')
p$add_argument('--outdir',             type="character",     help='File')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$diff_results_dir <- file.path(io$basedir,"results/atac/archR/differential/metacells/celltype/PeakMatrix")
# args$outdir <- file.path(io$basedir,"results/atac/archR/differential/metacells/celltype/PeakMatrix/parsed")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

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

fwrite(rbindlist(diff_results_list), file.path(args$outdir,"diff_results.txt.gz"), sep="\t", quote=F, na="NA")

##########
## TEST ##
##########

tmp <- fread("/bi/group/reik/ricard/data/gastrulation_multiome_10x/test/results/atac/archR/differential/pseudobulk/celltype/PeakMatrix/parsed/diff_results.txt.gz") %>%
  .[abs(logFC)>=2 & padj_fdr<=0.01 & (mean_groupA>=2.5 | mean_groupB>=2.5)] %>% .[,.N,by="feature"]

sum(tmp$N>=2)