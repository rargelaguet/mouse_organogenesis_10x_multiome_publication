source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/rna/cell_metadata.txt.gz")
io$indir <- file.path(io$basedir,"results/rna/cell_cycle")
io$outfile <- file.path(io$basedir,"results/rna/cell_cycle/cell_metadata_cell_cycle.txt.gz")

# Options

#################################################
## Load cell cycle assignments for each sample ##
#################################################

cell_cycle.dt <- opts$samples %>% map(function(i) {
  fread(file.path(io$indir,sprintf("%s_cell_cycle.txt.gz",i))) %>% .[,sample:=i]
}) %>% rbindlist

table(cell_cycle.dt$sample)

stopifnot(unique(cell_cycle.dt$sample)%in%opts$samples)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[sample%in%opts$samples]

##############################
## Merge with cell metadata ##
##############################

stopifnot(cell_cycle.dt$cell%in%cell_metadata.dt$cell)

cell_metadata_cell_cycle.dt <- cell_metadata.dt %>% merge(cell_cycle.dt, by=c("cell","sample"), all.x=TRUE)

# Rename columns
# Phase = cell_cycle_phase
# S.Score = cell_cycle_S_score
# G2M.Score = cell_cycle_G2M_score

##########
## Save ##
##########

fwrite(cell_metadata_cell_cycle.dt, io$outfile, sep="\t", quote=F)