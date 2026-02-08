  
source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"rna/results/cell_cycle/sample_metadata_with_cell_cycle_scores.txt")
io$outfile <- file.path(io$basedir,"rna/sample_metadata_standarised.txt")

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata)

sample_metadata.dt[,c("cell_line_id","distribution_bank"):=NULL]

#############################
## Update metadata columns ##
#############################

sample_metadata.dt %>% 
  setnames("sample","sample_name_original") %>%
  setnames("alias","sample")

# new cols
sample_metadata.dt %>%
  .[,biological_sex := "female"] %>%
  .[,tissue := "lung"] %>%
  .[,cell_line := "WI-38"] %>%
  .[,cell_type := "fibroblast"]


generic_cols <- c("sample", "age", "age_unit", "biological_sex", "batch", "disease", "ethnicity", "genotype", "tissue", "cell_line", "cell_type", "genetic_perturbation_target_gene", "genetic_perturbation_action", "chemical_perturbation_chemical", "chemical_perturbation_dose", "chemical_perturbation_dose_units and chemical_perturbation_time_duration")
specific_cols <- c("sample_name_original","condition", "time", "cell_cycle_score")

desired_cols <- c(generic_cols,specific_cols)
sample_metadata.dt <- sample_metadata.dt[,colnames(sample_metadata.dt)%in%desired_cols,with=F]

# Add missing columns with NA values
missing_cols <- setdiff(generic_cols, names(sample_metadata.dt))
if (length(missing_cols) > 0) {
  sample_metadata.dt[, (missing_cols) := NA_character_] # Use NA_character_ for character columns, or NA for mixed types
}

# Set the column order
setcolorder(sample_metadata.dt, desired_cols)

# Display the final data.table structure
print(sample_metadata.dt)
str(sample_metadata.dt)

##########
## Save ##
##########

fwrite(sample_metadata.dt, io$outfile, sep="\t", quote=F, na="NA")
