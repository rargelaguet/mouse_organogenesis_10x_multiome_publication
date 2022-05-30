here::i_am("atac/archR/differential/compare_differential_atac_pseudobulk_metacells_cells.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$matrix <- "PeakMatrix"  # "GeneScoreMatrix_TSS"
opts$group_variable <- "genotype"
opts$celltypes <- c("NMP") # only two groups

# I/O
io$basedir <- file.path(io$basedir,"test")
io$diff.pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/differential/pseudobulk/%s/%s",opts$group_variable, opts$matrix))
io$diff.cells <- file.path(io$basedir,sprintf("results/atac/archR/differential/cells/%s/%s",opts$group_variable, opts$matrix))
io$diff.metacells <- file.path(io$basedir,sprintf("results/atac/archR/differential/metacells/%s/%s",opts$group_variable, opts$matrix))
io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/comparison/%s/%s",opts$group_variable, opts$matrix)); dir.create(io$outdir, showWarnings = F, recursive = T)

# /Users/argelagr/data/gastrulation_multiome_10x/test/results/atac/archR/differential/metacells/genotype/PeakMatrix
# /Users/argelagr/data/gastrulation_multiome_10x/test/results/atac/archR/differential/pseudobulk/genotype/PeakMatrix
# /Users/argelagr/data/gastrulation_multiome_10x/test/results/atac/archR/differential/cells/genotype/PeakMatrix

##########################################
## Load results at the pseudobulk level ##
##########################################

# i <- "Epiblast"; j <- "Primitive_Streak"
atac_diff_pseudobulk.dt <- opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_WT_vs_KO_pseudobulk.txt.gz", io$diff.pseudobulk,j)
  if (file.exists(file)) {
    fread(file, select=c(1,2)) %>% 
      setnames("idx","feature") %>%
      .[,celltype:=j] %>%
      .[,class:="pseudobulk"]
  }
}) %>% rbindlist
  

########################################
## Load results at the metacell level ##
########################################

atac_diff_metacells.dt <- opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_WT_vs_KO.txt.gz", io$diff.metacells,j)
  if (file.exists(file)) {
    fread(file, select=c(1,2)) %>% 
      setnames("logFC","diff") %>%
      .[,celltype:=j] %>%
      .[,class:="metacells"]
  }
}) %>% rbindlist

# temporary
# atac_diff_metacells.dt %>% .[is.na(diff),diff:=0]
  
####################################
## Load results at the cell level ##
####################################

# j <- "Epiblast"; i <- "Primitive_Streak"
atac_diff_cells.dt <- opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_WT_vs_KO.txt.gz", io$diff.cells,j)
  if (file.exists(file)) {
    fread(file, select=c(1,2)) %>% 
      setnames(c("feature","diff")) %>%
      .[,celltype:=j] %>%
      .[,class:="cells"]
  }
}) %>% rbindlist

###################
## Sanity checks ##
###################

# all(sort(unique(atac_diff_cells.dt$peak))==sort(unique(atac_diff_metacells.dt$peak)))
# all(sort(unique(atac_diff_cells.dt$peak))==sort(unique(atac_diff_pseudobulk.dt$peak)))
# mean(is.na(atac_diff_metacells.dt$diff))
# mean(is.na(atac_diff_cells.dt$diff))
# mean(is.na(atac_diff_pseudobulk.dt$diff))

###########
## Merge ##
###########

stopifnot(colnames(atac_diff_pseudobulk.dt)==colnames(atac_diff_cells.dt))
stopifnot(colnames(atac_diff_pseudobulk.dt)==colnames(atac_diff_metacells.dt))

atac_diff.dt <- rbindlist(list(atac_diff_cells.dt, atac_diff_metacells.dt, atac_diff_pseudobulk.dt)) %>%
  dcast(feature+celltype~class, value.var="diff")

##########
## Plot ##
##########

to.plot <- atac_diff.dt %>% 
  .[abs(metacells)>=0.25] %>%
  .[sample.int(.N,size=1e4)]

ggscatter(to.plot, x="cells", y="pseudobulk", size=0.5, add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="Differential acc. (cells)", y="Differential acc. (pseudobulk)")

ggscatter(to.plot, x="cells", y="metacells", size=0.5, add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="Differential acc. (cells)", y="Differential acc. (metacells)")

ggscatter(to.plot, x="pseudobulk", y="metacells", size=0.5, add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="Differential acc. (pseudobulk)", y="Differential acc. (metacells)")
