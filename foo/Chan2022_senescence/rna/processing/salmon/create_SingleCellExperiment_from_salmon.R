source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$input <- file.path(io$basedir,"rna/processed/salmon/salmon.merged.gene_counts.tsv.gz")
io$outfile <- file.path(io$basedir,"rna/processed/salmon/salmon_genes_SingleCellExperiment.rds")

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata)

###############
## Load data ##
###############

counts.dt <- fread(io$input) 

colnames(counts.dt) <- gsub("_Homo_sapiens_RNA.Seq","",colnames(counts.dt))

counts.dt <- counts.dt[!duplicated(gene_name)]

gene_metadata.dt <- counts.dt[,c("gene_id","gene_name")]
counts.mtx <- counts.dt[,-c("gene_name")] %>% matrix.please %>% round
colnames(counts.mtx) <- colnames(counts.mtx)# %>% gsub("\\.","-",.)

sample_metadata.dt$sample[!sample_metadata.dt$sample%in%colnames(counts.mtx)]
counts.mtx <- counts.mtx[,sample_metadata.dt$sample]

#########################################
## Create SummarizedExperiment for TEs ##
#########################################

coldata_to_se <- DataFrame(sample_metadata.dt %>% tibble::column_to_rownames("sample"))
stopifnot(colnames(counts.mtx)==rownames(coldata_to_se))

rowdata_to_se <- DataFrame(gene_metadata.dt %>% tibble::column_to_rownames("gene_id"))
stopifnot(rownames(counts.mtx)==rownames(rowdata_to_se))

sce <- SingleCellExperiment(
  assays = list(counts=counts.mtx), 
  rowData = rowdata_to_se,
  colData = coldata_to_se
)

##########
## Save ##
##########

saveRDS(sce, io$outfile)
