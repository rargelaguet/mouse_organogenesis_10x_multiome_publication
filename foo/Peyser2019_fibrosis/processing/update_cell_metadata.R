
source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# io$metadata <- file.path(io$basedir,"processed/cell_metadata.txt.gz")
# io$outfile <- file.path(io$basedir,"processed/cell_metadata.txt.gz")

io$metadata <- file.path(io$basedir,"results/metacells/per_condition/metacells_metadata.txt.gz")
io$outfile <- file.path(io$basedir,"results/metacells/per_condition/metacells_metadata.txt.gz")

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(io$metadata)

#####################
## Update metadata ##
#####################

# cell_metadata.dt[,condition_original:=condition]
# cell_metadata.dt[condition=="Saline Control",condition:="Saline_Control"]
# cell_metadata.dt[condition=="Bleomycin Fibrosis Model",condition:="Bleomycin_day11"]
# cell_metadata.dt[condition=="Bleomycin Fibrosis plus Nintedanib",condition:="Bleomycin_nintedanib_day11"]
# unique(cell_metadata.dt$condition)

cell_metadata.dt[,cell_type:=gsub(" ","_",cell_type)]
unique(cell_metadata.dt$cell_type)

##########
## Save ##
##########

fwrite(cell_metadata.dt, io$outfile, sep="\t", quote=F, na="NA")
