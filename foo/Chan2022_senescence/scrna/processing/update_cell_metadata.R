
source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/cell_cycle/cell_metadata_cell_cycle.txt.gz")
io$outfile <- file.path(io$basedir,"results/cell_cycle/cell_metadata_cell_cycle.txt.gz")

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(io$metadata)

#####################
## Update metadata ##
#####################

cell_metadata.dt %>% .[,day:=opts$sample2day[sample]]

cell_metadata.dt[,c("alias","day")] %>% unique

##########
## Save ##
##########

fwrite(cell_metadata.dt, io$outfile, sep="\t", quote=F)
