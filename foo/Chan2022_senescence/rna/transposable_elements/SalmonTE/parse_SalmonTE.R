source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$inputdir <- file.path(io$basedir,"SalmonTE")
io$outdir <- file.path(io$basedir,"processed/salmonTE"); dir.create(io$outdir, showWarnings=F)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata)

opts$samples <- sample_metadata.dt$sample

###############
## Load data ##
###############

stopifnot(dir.exists(file.path(io$inputdir,opts$samples)))

TE_expression_list <- list();
          
# i <- "SRR9126653"
for (i in opts$samples) {
  TE_expression_list[[i]] <- fread(file.path(io$inputdir,sprintf("%s/quant.sf",i))) %>%
    .[,c("Name","NumReads")] %>%
    setnames(c("feature",i)) %>% 
    matrix.please
}

###########
## Merge ##
###########

TE_expression.mtx <- do.call("cbind",TE_expression_list)

########################
## Create TE metadata ##
########################

TE_metadata.dt <- fread(file.path(io$basedir,"SalmonTE/clades.csv")) %>%
  .[,class:=NULL] %>% setnames(c("name","class"))

features.to.use <- intersect(rownames(TE_expression.mtx),TE_metadata.dt$name)

# filter
TE_metadata.dt <- TE_metadata.dt[name%in%features.to.use]
TE_expression.mtx <- TE_expression.mtx[TE_metadata.dt$name,]

#########################################
## Create SummarizedExperiment for TEs ##
#########################################

coldata_to_se <- DataFrame(sample_metadata.dt[sample%in%opts$samples] %>% setkey(sample) %>% .[colnames(TE_expression.mtx)] %>% tibble::column_to_rownames("sample"))
stopifnot(colnames(TE_expression.mtx)==rownames(coldata_to_se))

rowdata_to_se <- DataFrame(TE_metadata.dt %>% tibble::column_to_rownames("name"))
stopifnot(rownames(TE_expression.mtx)==rownames(rowdata_to_se))

te.se <- SingleCellExperiment(
  assays = list(counts=TE_expression.mtx), 
  rowData = rowdata_to_se,
  colData = coldata_to_se
)

##########
## Save ##
##########

saveRDS(te.se, file.path(io$outdir,"SalmonTE_SingleCellExperiment.rds"))
