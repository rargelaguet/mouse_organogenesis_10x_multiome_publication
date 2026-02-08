source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
io$counts <- file.path(io$basedir,"processed/atac/coverage_peaks/atac_coverage_peaks.count.gz")
io$outfile <- file.path(io$basedir,"processed/atac/SingleCellExperiment.rds")

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata)

samples <- as.character(read.table("/share/rargelaguet/data/Chan2022_senescence/processed/atac/coverage_peaks/sample_names.txt")[1,])
stopifnot(samples%in%sample_metadata.dt$sample)

#################
## Load counts ##
#################

counts.dt <- fread(io$counts) %>% 
  setnames(c("chr","start","end",samples)) %>%
  .[,chr:=paste0("chr",chr)]

#############################
## Define feature metadata ##
#############################

feature_metadata.dt <- counts.dt[,c("chr","start","end")] %>% 
  .[,feature:=sprintf("%s:%s-%s",chr,start,end)]

#########################
## Define count matrix ##
#########################

counts.mtx <- counts.dt %>% 
  .[,c("chr","start","end"):=NULL] %>% 
  matrix.please

rownames(counts.mtx) <- feature_metadata.dt$feature

#################################
## Create SingleCellExperiment ##
#################################

coldata_to_sce <- DataFrame(sample_metadata.dt[sample%] %>% setkey(sample) %>% .[colnames(counts.mtx)] %>% tibble::column_to_rownames("sample"))
stopifnot(colnames(counts.mtx)==rownames(coldata_to_sce))

rowdata_to_sce <- DataFrame(feature_metadata.dt %>% tibble::column_to_rownames("feature"))
stopifnot(rownames(counts.mtx)==rownames(rowdata_to_sce))

sce <- SingleCellExperiment(
  assays = list(counts=counts.mtx), 
  rowData = rowdata_to_sce,
  colData = coldata_to_sce
)

##########
## Save ##
##########

saveRDS(sce, io$outfile)
