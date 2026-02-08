  
source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"cell_metadata.txt.gz")
io$outfile <- file.path(io$basedir,"cell_metadata_standarised.txt.gz")

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(io$metadata) %>%
  .[,c("cell", "alias", "group", "genetic_background", "ncounts", "nfeatures", "mitochondrial_percent_RNA", "cell_cycle_phase")] %>%
  setnames(c("cell", "sample", "group", "genetic_background", "number_reads", "number_features", "mito_percentage", "cell_cycle_phase"))

cell_metadata.dt[grep("PDL_50",sample),c("cell_cycle_phase"):=list("G1")]

table(cell_metadata.dt$group)

#############################
## Update metadata columns ##
#############################

# new cols
cell_metadata.dt %>%
  .[,cell_type := "fibroblast"] %>%
  .[,biological_sex := "female"] %>%
  .[,tissue := "lung"] %>%
  .[,cell_line := "WI-38"]


# Set the column order
desired_cols <- c("cell", "sample", "group", "genetic_background", "cell_type", "biological_sex", "tissue", "cell_line", "cell_cycle_phase", "number_reads", "number_features", "mito_percentage")
setcolorder(cell_metadata.dt, desired_cols)

# Display the final data.table structure
print(cell_metadata.dt)
str(cell_metadata.dt)

##########
## Save ##
##########

fwrite(cell_metadata.dt, io$outfile, sep="\t", quote=T, na="NA")
