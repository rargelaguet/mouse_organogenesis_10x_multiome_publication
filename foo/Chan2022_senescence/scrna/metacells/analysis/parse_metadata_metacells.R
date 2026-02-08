source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata.cells <- file.path(io$basedir,"cell_metadata.txt.gz")
io$metadata.metacells <- file.path(io$basedir,"results/rna/metacells/scTE/metacells_metadata.txt.gz")
io$cell2metacell_assignment <- file.path(io$basedir,"results/rna/metacells/scTE/cell2metacell_assignment.txt.gz")
# io$outdir <- paste0(io$basedir,"/results/metacells")

# Options

############################
## Load metacell metadata ##
############################

metacell_metadata.dt <- fread(io$metadata.metacells) 

cell_metadata.dt <- fread(io$metadata.cells)

cell2metacell.dt <- fread(io$cell2metacell_assignment)
  
####################################################
## Add cell type assignments (mapping to Liu2020) ##
####################################################

# io$mapping <- file.path(io$basedir,"results/mapping_liu2020/cell_metadata_after_mapping_liu2020.txt.gz")
# 
# mapping_liu2020.dt <- fread(io$mapping) %>%
#   .[,c("cell","cell_liu2020","alias_liu2020","clustering_liu2020","celltype_liu2020")] %>%
#   setnames("clustering_liu2020","cluster_liu2020")
# 
# tmp <- cell2metacell.dt %>%
#   merge(mapping_liu2020.dt[,c("cell","cluster_liu2020")],by="cell") %>%
#   .[,.N,by=c("metacell","cluster_liu2020")] %>%
#   .[,.SD[which.max(N)],by="metacell"] %>% .[,N:=NULL]
# metacell_metadata.dt <- metacell_metadata.dt %>% merge(tmp,by="metacell")
# 
# opts$cluster2celltype <- c(
#   "fm1" = "Fibroblast_adult", 
#   "fm2" = "Fibroblast_adult", 
#   "fm3" = "Fibroblast_reprogrammed", 
#   "fm4" = "Fibroblast_reprogrammed", 
#   "fm5" = "Fibroblast_reprogrammed", 
#   "mix" = "Mixed_state", 
#   "re1" = "Refractory", 
#   "re2" = "Refractory", 
#   "re3" = "Refractory", 
#   "re4" = "Refractory", 
#   "re5" = "Refractory", 
#   "re6" = "Refractory", 
#   "nic" = "nic", 
#   "pr1" = "Primed_iPSCs", 
#   "pr2" = "Primed_iPSCs", 
#   "pr3" = "Primed_iPSCs", 
#   "pr4" = "Primed_iPSCs", 
#   "nr1" = "Naive_iPSCs", 
#   "nr2" = "Naive_iPSCs", 
#   "nr3" = "Naive_iPSCs", 
#   "nr4" = "Naive_iPSCs"
# )
# metacell_metadata.dt[,celltype_liu2020:=str_replace_all(cluster_liu2020,opts$cluster2celltype)]

##############################
## Update metacell metadata ##
##############################

tmp <- cell2metacell.dt %>% 
  merge(cell_metadata.dt[,c("cell","cluster_liu2020","celltype_liu2020")],by="cell") %>%
  .[,.N,by=c("metacell","cluster_liu2020","celltype_liu2020")] %>% 
  .[,.SD[which.max(N)],by="metacell"] %>% .[,N:=NULL]

metacell_metadata.dt <- metacell_metadata.dt %>% merge(tmp,by="metacell")

##########
## save ##
##########

fwrite(metacell_metadata.dt, io$metadata.metacells, sep="\t", na="NA", quote=F)
