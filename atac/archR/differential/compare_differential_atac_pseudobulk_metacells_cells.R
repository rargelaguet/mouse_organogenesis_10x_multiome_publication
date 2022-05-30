here::i_am("atac/archR/differential/compare_differential_atac_pseudobulk_metacells_cells.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$matrix <- "PeakMatrix"
# opts$matrix <- "GeneScoreMatrix_TSS"
opts$group_variable <- "celltype"
opts$groups <- c("Erythroid3", "Neural_crest") # only two groups

# I/O
io$basedir <- file.path(io$basedir,"test")
io$diff.pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/differential/pseudobulk/%s/%s",opts$group_variable, opts$matrix))
io$diff.cells <- file.path(io$basedir,sprintf("results/atac/archR/differential/cells/%s/%s",opts$group_variable, opts$matrix))
io$diff.metacells <- file.path(io$basedir,sprintf("results/atac/archR/differential/metacells/%s/%s",opts$group_variable, opts$matrix))
io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/comparison/%s/%s",opts$group_variable, opts$matrix)); dir.create(io$outdir, showWarnings = F, recursive = T)

##########################################
## Load results at the pseudobulk level ##
##########################################

# i <- "Epiblast"; j <- "Primitive_Streak"
atac_diff_pseudobulk.dt <- opts$groups %>% map(function(i) { opts$groups %>% map(function(j) {
  file <- file.path(io$diff.pseudobulk,sprintf("%s_vs_%s_pseudobulk.txt.gz",i,j))
  if (file.exists(file)) {
    fread(file, select = c(1,2)) %>% 
      setnames(c("peak","diff")) %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      .[,class:="pseudobulk"] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist

########################################
## Load results at the metacell level ##
########################################

atac_diff_metacells.dt <- opts$groups %>% map(function(i) { opts$groups %>% map(function(j) {
  file <- file.path(io$diff.metacells,sprintf("%s_%s_vs_%s.txt.gz",opts$matrix,i,j))
  if (file.exists(file)) {
    fread(file, select = c(1,2)) %>% 
      setnames(c("peak","diff")) %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      .[,class:="metacells"] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist

# temporary
atac_diff_metacells.dt %>% .[is.na(diff),diff:=0]
  
####################################
## Load results at the cell level ##
####################################

# j <- "Epiblast"; i <- "Primitive_Streak"
atac_diff_cells.dt <- opts$groups %>% map(function(i) { opts$groups %>% map(function(j) {
  file <- file.path(io$diff.cells,sprintf("%s_%s_vs_%s.txt.gz",opts$matrix,i,j))
  if (file.exists(file)) {
    fread(file, select=c(1,2)) %>% 
      setnames(c("peak","diff")) %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      .[,class:="cells"] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist

# ad hoc
# if ("name"%in%colnames(atac_diff_cells.dt)) {
#   atac_diff_cells.dt[,idx:=NULL] %>% setnames("name","idx")
# }

###################
## Sanity checks ##
###################

# all(sort(unique(atac_diff_cells.dt$peak))==sort(unique(atac_diff_metacells.dt$peak)))
# all(sort(unique(atac_diff_cells.dt$peak))==sort(unique(atac_diff_pseudobulk.dt$peak)))
# mean(is.na(atac_diff_metacells.dt$diff))
# mean(is.na(atac_diff_cells.dt$diff))
# mean(is.na(atac_diff_pseudobulk.dt$diff))

# tmp <- rbindlist(list(atac_diff_cells.dt, atac_diff_metacells.dt, atac_diff_pseudobulk.dt))
# tmp[,.N,by=c("class","celltypeA","celltypeB")]

if (!atac_diff_cells.dt$celltypeA[1]==opts$groups[1]) {
  atac_diff_cells.dt %>% .[,diff:=-diff] %>% .[,c("celltypeA","celltypeB"):=list(celltypeB,celltypeA)]
}
if (!atac_diff_metacells.dt$celltypeA[1]==opts$groups[1]) {
  atac_diff_metacells.dt %>% .[,diff:=-diff] %>% .[,c("celltypeA","celltypeB"):=list(celltypeB,celltypeA)]
}
if (!atac_diff_pseudobulk.dt$celltypeA[1]==opts$groups[1]) {
  atac_diff_pseudobulk.dt %>% .[,diff:=-diff] %>% .[,c("celltypeA","celltypeB"):=list(celltypeB,celltypeA)]
}

###########
## Merge ##
###########

atac_diff.dt <- rbindlist(list(atac_diff_cells.dt, atac_diff_metacells.dt, atac_diff_pseudobulk.dt)) %>%
  dcast(peak+celltypeA+celltypeB~class, value.var="diff")

##########
## Plot ##
##########

to.plot <- atac_diff.dt %>% .[sample.int(.N,size=1e4)]

ggscatter(to.plot, x="cells", y="pseudobulk", size=0.5, add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="Differential acc. (cells)", y="Differential acc. (pseudobulk)")

ggscatter(to.plot, x="cells", y="metacells", size=0.5, add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="Differential acc. (cells)", y="Differential acc. (metacells)")

ggscatter(to.plot, x="pseudobulk", y="metacells", size=0.5, add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="Differential acc. (pseudobulk)", y="Differential acc. (metacells)")
