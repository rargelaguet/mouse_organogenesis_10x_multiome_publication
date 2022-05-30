here::i_am("rna/differential/other/extract_TFs_diff.R")

# Load default settings
source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--TFs',             type="character",     help='Cell metadata file')
p$add_argument('--diff_results',   type="character",     help='File')
p$add_argument('--outfile',             type="character",     help='File')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$TFs <- "/Users/argelagr/data/mm10_regulation/TFs/TFs.txt"
# args$diff_results <- file.path(io$basedir,"results/rna/differential/metacells/celltype/parsed/diff_expr_results.txt.gz")
# args$outfile <- file.path(io$basedir,"results/rna/differential/metacells/celltype/parsed/diff_expr_results_tfs.txt.gz")
## END TEST ##

dir.create(dirname(args$outfile), showWarnings = F)

##############
## Load TFs ##
##############

# TFs <- fread(args$TFs)[["gene"]]
TFs <- fread(args$TFs)[[1]] %>% str_to_title

################################################
## Load differential expression and fetch TFs ##
################################################

diff_results.dt <- fread(args$diff_results) %>% 
  .[gene%in%TFs] %>% .[,gene:=toupper(gene)]

# diff_tf.dt <- opts$celltypes %>% map(function(i) {
#   opts$celltypes %>% map(function(j) {
#     if (i!=j) {
#       file <- file.path(args$diff_results_dir,sprintf("%s_vs_%s.txt.gz",i,j))
#       if (file.exists(file)) {
#         fread(file, select=c(1,2,4,6,7,8,9)) %>% 
#           setnames(c("gene", "logFC", "padj_fdr", "groupA_N", "groupB_N", "groupA_detection_rate", "groupB_detection_rate")) %>%
#           .[gene%in%TFs] %>% .[,c("celltypeA","celltypeB"):=list(i,j)]
#       }
#     }
#   }) %>% rbindlist
# }) %>% rbindlist %>% .[,gene:=toupper(gene)]

print(sprintf("Number of TFs in the differential expression results: %s",length(unique(diff_results.dt$gene))))

##########
## Save ##
##########

fwrite(diff_results.dt, args$outfile, sep="\t", quote=F, na="NA")

