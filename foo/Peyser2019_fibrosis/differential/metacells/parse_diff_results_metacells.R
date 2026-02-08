source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# Options
opts$group_variable <- "condition"; opts$groups <- opts$conditions

# I/O
io$diff_results_dir <- file.path(io$basedir,sprintf("results/differential/metacells/%s",opts$group_variable))
io$outdir <- file.path(io$basedir,sprintf("results/differential/metacells/%s/parsed",opts$group_variable)); dir.create(io$outdir, showWarnings = F, recursive=T)

#########################
## Load and parse data ##
#########################

diff_results_list <- list()

celltypes.to.use <- c("Type_II_alveolar", "dendritic_cells", "monocytes", "Type_I_alveolar", "neutrophils", "mesothelium", "endothelial", "macrophages", "B_cells", "fibroblasts")

# i <- opts$groups[1]; j <- opts$groups[2]
for (i in celltypes.to.use) {
  
  # Control vs Bleomycin
  groupA <- "Saline_Control"; groupB <- "Bleomycin_day11"
  file <- file.path(io$diff_results_dir,sprintf("%s_%s_vs_%s.txt.gz",i,groupA,groupB))
  if (file.exists(file)) {
    diff_results_list[[sprintf("%s_%s_vs_%s",i,groupA,groupB)]] <- fread(file) %>% .[,c("celltype","groupA","groupB"):=list(i,groupA,groupB)]
  } else {
    print(sprintf("%s not found...",file))
  }
  
  # Control vs Bleomycin+Nintedanib
  groupA <- "Saline_Control"; groupB <- "Bleomycin_nintedanib_day11"
  file <- file.path(io$diff_results_dir,sprintf("%s_%s_vs_%s.txt.gz",i,groupA,groupB))
  if (file.exists(file)) {
    diff_results_list[[sprintf("%s_%s_vs_%s",i,groupA,groupB)]] <- fread(file) %>% .[,c("celltype","groupA","groupB"):=list(i,groupA,groupB)]
  } else {
    print(sprintf("%s not found...",file))
  }
  
  # Bleomycin vs Bleomycin+Nintedanib
  groupA <- "Bleomycin_day11"; groupB <- "Bleomycin_nintedanib_day11"
  file <- file.path(io$diff_results_dir,sprintf("%s_%s_vs_%s.txt.gz",i,groupA,groupB))
  if (file.exists(file)) {
    diff_results_list[[sprintf("%s_%s_vs_%s",i,groupA,groupB)]] <- fread(file) %>% .[,c("celltype","groupA","groupB"):=list(i,groupA,groupB)]
  } else {
    print(sprintf("%s not found...",file))
  }
  
}

##########
## Save ##
##########

fwrite(rbindlist(diff_results_list), file.path(io$outdir,"diff_expr_results.txt.gz"), sep="\t", quote=F, na="NA")

