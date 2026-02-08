#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
io$se <- file.path(io$basedir,"processed/atac/jack/SummarizedExperiment_ATAC_Jack.rds")
io$outfile <- paste0(io$basedir,"/processed/atac/jack/SingleCellExperiment.rds")

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata)

sample_metadata.dt$alias <- opts$sample2alias[sample_metadata.dt$sample]

# sample_metadata.dt %>% setnames("TreatmentGroup","Group")

###############
## Load data ##
###############

se <- readRDS(io$se)[,sample_metadata.dt$sample]

# Add sample metadata as colData
colData(se) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

# Rename peaks
tmp <- as.data.table(rowRanges(se)) %>% .[,c("seqnames","start","end","name")] %>% .[,feature:=sprintf("%s:%s-%s",seqnames,start,end)]
rownames(se) <- tmp$feature

rowRanges(se)
# rownames(se) <- rowData(se)$gene_name

# Remove duplicated gene names
se <- se[!duplicated(rownames(se))]

assays(se)["normCounts"] <- NULL
assays(se)["normCounts_CCrm"] <- NULL
assayNames(se) <- c("counts")

#####################################
## Convert to SingleCellExperiment ##
#####################################

sce <- SingleCellExperiment(assays = assays(se), colData = colData(se), rowData = rowData(se))
sce

##########
## Save ##
##########

# fwrite(sample_metadata.dt, "/Users/rargelaguet/data/Chan2022_senescence/processed/atac/sample_metadata.txt")
saveRDS(sce, io$outfile)
