  
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


#############################
## Update metadata columns ##
#############################

sample_metadata.dt %>% 
  setnames("sample","sample name") %>%
  setnames("cell_cycle_score","G2M Score")
  
sample_metadata.dt %>%
  .[,`cell population doubling`:=as.numeric(NA)] %>%
  .[time=="TP1",`cell population doubling`:=20] %>%
  .[time=="TP2",`cell population doubling`:=25] %>%
  .[time=="TP3",`cell population doubling`:=28] %>%
  .[time=="TP4",`cell population doubling`:=33] %>%
  .[time=="TP5",`cell population doubling`:=37] %>%
  .[time=="TP6",`cell population doubling`:=46] %>%
  .[time=="TP7",`cell population doubling`:=45] %>%
  .[time=="TP8",`cell population doubling`:=50] %>%
  .[time=="TP9",`cell population doubling`:=52] %>%
  .[time=="TP10",`cell population doubling`:=53]
  
sample_metadata.dt[,`other perturbation`:=condition]

sample_metadata.dt[,`other perturbation dose`:=NA_character_]
sample_metadata.dt[,`other perturbation dose units`:=NA_character_]
sample_metadata.dt[,`other perturbation time duration`:=time]

  
sample_metadata.dt %>%
  .[,`experiment accession` := "geo:GSE175533"] %>%
  .[,`experiment name` := "chan2022_senescence"] %>%
  .[,`experiment title` := "Senescence in WI-38 fibroblasts"] %>%
  .[,`contact name` := "Ricard Argelaguet"] %>%
  .[,`contact email` := "rargelaguet@altoslabs.com"] %>%
  .[,`donor accession` := "coriell:AG06814-N"] %>%
  .[,`donor name` := "AG06814-N"] %>%
  .[,`donor genotype` := NA_character_] %>%
  .[,`donor species` := "Homo sapiens"] %>%
  .[,`donor biological sex` := "female"] %>%
  .[,`donor age` := 16] %>%
  .[,`donor age range` := NA] %>%
  .[,`donor age units` := "fetal week"] %>%
  .[,`donor developmental stage` := "fetal stage"] %>%
  .[,`human ethnicity` := "European"] %>%
  .[,`donor disease` := "normal"] %>%
  .[,`specimen accession` := "coriell:AG06814-N"] %>%
  .[,`specimen name` := "AG06814-N"] %>%
  .[,`specimen type` := "cell line"] %>%
  .[,`specimen tissue` := "lung primordium"] %>%
  # .[,`sample prep name` := NA_character_] %>%
  # .[,`batch` := NA_character_] %>%
  .[,`cell category` := "primary"] %>%  # CHANGE HTERT TO IMMORTALSIED
  .[,`cell line` := "WI-38"] %>%
  .[,`cell type` := "lung fibroblast"] %>%
  .[,`culture media` := "DMEM, 10% dialyzed FBS"] %>%
  .[,`biological replicate` := strsplit(alias,"rep") %>% map_chr(2)]




generic_cols <- c(
  "sample name",
  "experiment name",
  "experiment title",
  "donor accession",
  "donor name",
  "biological replicate",
  "donor genotype",
  "donor species",
  "donor biological sex",
  "donor age",
  "donor age range",
  "donor age units",
  "human ethnicity",
  "donor disease",
  # "sample prep name",
  "batch",
  "specimen tissue",
  "primary cell selection markers",
  "cell line",
  "cell type",
  "genetic perturbation name",
  "genetic perturbation target gene",
  "genetic perturbation action",
  "chemical perturbation name",
  "chemical perturbation chemical",
  "chemical perturbation chemical type",
  "chemical perturbation dose",
  "chemical perturbation dose units",
  "chemical perturbation time duration",
  "chemical perturbation start time",
  "chemical perturbation end time",
  "chemical perturbation time units",
  "other perturbation",
  "other perturbation dose",
  "other perturbation dose units",
  "other perturbation time duration",
  "G2M Score"
)

specific_cols <- c()
desired_cols <- c(specific_cols,generic_cols)
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
