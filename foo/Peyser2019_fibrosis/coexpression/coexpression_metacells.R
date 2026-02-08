
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))

set.seed(42)

#####################
## Define settings ##
#####################

# I/O
io$sce <- file.path(io$basedir,"results/rna/metacells/scTE/SingleCellExperiment_metacells.rds")
io$metadata <- file.path(io$basedir,"results/rna/metacells/scTE/metacells_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/coexpression"); dir.create(io$outdir, showWarnings = F)

# Options
opts$groups <- c(
  # Baseline
  # "25F_day0",
  # "28F_day0",
  # "52F_day0",
  # "54F_day0",
  "63F_day0",
  # "70M_day0",
  # "FSFA_day0",
  "FSFB_day0",
  # 63F OSKM
  "63F_OSKM-T3G_day2",
  "63F_OSKM-T3G_day7",
  "63F_OSKM-T3G_day13",
  "63F_OSKM-T3G_day27_unsorted",
  # BJ_CHAMP OSKM
  "BJ_OSKM-CHAMP_day0",
  "BJ_OSKM-CHAMP_day2",
  "BJ_OSKM-CHAMP_day7",
  "BJ_OSKM-CHAMP_day13",
  "BJ_OSKM-CHAMP_day27_FTR",
  # FSF OSKM
  "FSF_OSKM-T3G_day2",
  "FSF_OSKM-T3G_day7",
  "FSF_OSKM-T3G_day13",
  "FSF_OSKM-T3G_day27_FTR",
  "FSF_OSKM-T3G_day27_TR",
  "FSF_OSKM-T3G_day27_unsorted"
)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[group%in%opts$groups]

table(cell_metadata.dt$group)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$metacell, normalise = T)

# Add sample metadata as colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

###########################
## Coexpression analysis ##
###########################

tmp <- apply(logcounts(sce),1,var)
genes.to.use <- names(which(tmp>0.25))

mat <- logcounts(sce)[genes.to.use,]
mat = mat - rowMeans(mat)
mat = mat / sqrt(rowSums(mat^2))

cor.mtx = tcrossprod(mat) %>% round(2)

##########
## Save ##
##########

saveRDS(cor.mtx, file.path(io$outdir,"correlation_matrix.rds"))

