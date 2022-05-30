# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

source(here::here("rna/mapping/analysis/plot_utils.R"))

#####################
## Define settings ##
#####################

## I/O
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacells_metadata.txt.gz")
io$metacell_sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")
io$trajectory <- file.path(io$basedir,"results/rna/trajectories/nmp/nmp_trajectory.txt.gz")
io$atlas_trajectory <- file.path(io$atlas.basedir,"results/trajectories/nmp_somitic_spinal/nmp_trajectory.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/pdf"); dir.create(io$outdir, showWarnings = F)

# Dot size
opts$size.mapped <- 1
opts$size.nomapped <- 0.1

# Dot transparency
opts$alpha.mapped <- 0.85
opts$alpha.nomapped <- 0.35

###################
## Load metadata ##
###################

metacell_metadata.dt <- fread(io$metacell_metadata)
# sample_metadata.dt <- fread(io$metadata)

###############################
## Load SingleCellExperiment ##
###############################

sce <- readRDS(io$metacell_sce)

#####################
## Load trajectory ##
#####################

# trajectory.dt <- fread(io$atlas_trajectory) %>% setnames(c("cell","V1","V2"))
trajectory.dt <- fread(io$trajectory) %>% setnames(c("cell","V1","V2"))

#################################################
## Plot mapping of metacells to the trajectory ##
#################################################

to.plot <- trajectory.dt %>% copy %>%
  # .[,index:=match(cell, metacell_metadata.dt$closest.cell)] %>% 
  .[,index:=match(cell, metacell_metadata.dt$metacell)] %>% 
  .[,mapped:=as.factor(!is.na(index))] %>% 
  .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas","Metacell"))] %>%
  setorder(mapped) 

p <- plot.dimred(to.plot, query.label = "Metacell", atlas.label = "Atlas")

pdf(file.path(io$outdir,"trajectory_highlight_metacells.pdf"), width=8, height=6.5)
print(p)
dev.off()
