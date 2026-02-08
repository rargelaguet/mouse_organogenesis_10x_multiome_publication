
source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"cell_metadata.txt.gz")
io$outfile <- file.path(io$basedir,"sample_metadata.txt")

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(io$metadata)

############################
## Create sample metadata ##
############################

cols <- c("sample", "alias", "group", "day", "genetic_background")

sample_metadata.dt <- cell_metadata.dt[,cols,with=F] %>% unique

###################################
## Update cell cycle information ##
###################################

cell_cycle_stats.dt <- cell_metadata.dt[, .(
  g1_fraction = mean(cell_cycle_phase == "G1", na.rm = TRUE) %>% round(2),
  g2m_fraction = mean(cell_cycle_phase == "G2M", na.rm = TRUE) %>% round(2),
  s_fraction = mean(cell_cycle_phase == "S", na.rm = TRUE) %>% round(2)
), by = .(sample)]

sample_metadata.dt <- merge(sample_metadata.dt, cell_cycle_stats.dt, by = "sample", all.x = TRUE)

##########
## Save ##
##########

fwrite(sample_metadata.dt, io$outfile, sep="\t", quote=F)
