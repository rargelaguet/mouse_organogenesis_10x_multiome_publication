source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/cellranger/cell_metadata.txt.gz")
io$indir <- file.path(io$basedir,"results/cell_cycle")
io$outfile <- file.path(io$basedir,"results/cell_cycle/cell_metadata_cell_cycle.txt.gz")

# Options

#################################################
## Load cell cycle assignments for each sample ##
#################################################

cell_cycle.dt <- opts$samples %>% map(function(i) {
  fread(file.path(io$indir,sprintf("%s_cell_cycle.txt.gz",i))) %>% .[,sample:=i]
}) %>% rbindlist

# cell_cycle.dt[grep("PDL_50",sample),c("cell_cycle_s_score","cell_cycle_g2m_score","cell_cycle_phase"):=list(NA,NA,"G1")]
# table(cell_cycle.dt$sample)

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

setnames(cell_metadata_cell_cycle.dt, old = c("Phase", "S.Score", "G2M.Score"), new = c("cell_cycle_phase", "cell_cycle_S_score", "cell_cycle_G2M_score"))

##########
## Save ##
##########

fwrite(cell_metadata_cell_cycle.dt, io$outfile, sep="\t", quote=F)