here::i_am("atac/archR/chromvar_chip/metacells/differential/celltype/analysis/compare_differential_chromvar_pseudobulk_metacells.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$celltypes <- c("NMP","Epiblast","Gut")

# I/O
io$basedir <- file.path(io$basedir,"test")
io$diff.pseudobulk <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk/differential/celltypes/CISBP")
io$diff.metacells <- file.path(io$basedir,"results/atac/archR/chromvar_chip/metacells/differential/celltypes/CISBP")
# io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/comparison/%s/%s",opts$group_variable, opts$matrix)); dir.create(io$outdir, showWarnings = F, recursive = T)

##########################################
## Load results at the pseudobulk level ##
##########################################

# i <- "Epiblast"; j <- "Primitive_Streak"
chromvar_diff_pseudobulk.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$diff.pseudobulk,sprintf("chromVAR_%s_vs_%s.txt.gz",i,j))
  if (file.exists(file)) {
    fread(file, select=c(1,2)) %>% 
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>% 
      .[,class:="pseudobulk"] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist 


########################################
## Load results at the metacell level ##
########################################

chromvar_diff_metacells.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$diff.metacells,sprintf("chromVAR_%s_vs_%s.txt.gz",i,j))
  if (file.exists(file)) {
    fread(file, select=c(1,2)) %>% 
      .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>% 
      .[,class:="metacells"] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist 
  
###################
## Sanity checks ##
###################

# all(sort(unique(atac_diff_cells.dt$peak))==sort(unique(chromvar_diff_metacells.dt$peak)))
# all(sort(unique(atac_diff_cells.dt$peak))==sort(unique(chromvar_diff_pseudobulk.dt$peak)))
# mean(is.na(chromvar_diff_metacells.dt$diff))
# mean(is.na(atac_diff_cells.dt$diff))
# mean(is.na(chromvar_diff_pseudobulk.dt$diff))

###########
## Merge ##
###########

stopifnot(colnames(chromvar_diff_pseudobulk.dt)==colnames(chromvar_diff_metacells.dt))

chromvar_diff.dt <- rbindlist(list(chromvar_diff_metacells.dt, chromvar_diff_pseudobulk.dt)) %>%
  dcast(gene+celltypeA+celltypeB~class, value.var="diff")

##########
## Plot ##
##########

to.plot <- chromvar_diff.dt[celltypeA=="Epiblast" & celltypeB=="NMP"]

ggscatter(to.plot, x="pseudobulk", y="metacells", size=0.5, add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="Differential chromvar (pseudobulk)", y="Differential chromvar (metacells)")
