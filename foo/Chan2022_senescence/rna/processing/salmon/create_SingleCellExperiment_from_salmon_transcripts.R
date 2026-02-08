source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"processed/rna/sample_metadata.txt")
io$input <- file.path(io$basedir,"processed/rna/salmon/salmon.merged.transcript_tpm.tsv.gz")
io$outfile <- file.path(io$basedir,"processed/rna/salmon/salmon_transcripts_SingleCellExperiment.rds")

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata)

##########################
## Load gene metadata ##
##########################

gene_metadata.dt <- fread(io$gene_metadata)

###############
## Load data ##
###############

counts.dt <- fread(io$input) 

counts.dt <- counts.dt[!duplicated(gene_name)]

gene_metadata.dt <- counts.dt[,c("tx","gene_id")] %>% 
  setnames(c("transcript_id","ens_id")) %>% 
  merge(gene_metadata.dt[,c("ens_id","symbol")])

counts.mtx <- counts.dt[,-c("gene_id")] %>% matrix.please %>% round
colnames(counts.mtx) <- c(
  "Control_p8_rep1",
  "Control_p8_rep2",
  "Control_p8_rep3",
  "TG_p8_rep1",
  "TG_p8_rep2",
  "TG_p8_rep3",
  "Control_p11_rep1",
  "Control_p11_rep2",
  "Control_p11_rep3",
  "TG_p11_rep1",
  "TG_p11_rep2",
  "TG_p11_rep3",
  "Control_p14_rep1",
  "Control_p14_rep2",
  "Control_p14_rep3",
  "TG_p14_rep1",
  "TG_p14_rep2",
  "TG_p14_rep3",
  "Control_p19_rep1",
  "Control_p19_rep2",
  "Control_p19_rep3",
  "TG_p19_rep1",
  "TG_p19_rep2",
  "TG_p19_rep3",
  "Control_p22_rep1",
  "Control_p22_rep2",
  "Control_p22_rep3",
  "TG_p22_rep1",
  "TG_p22_rep2",
  "TG_p22_rep3",
  "Control_p24_rep1",
  "Control_p24_rep2",
  "Control_p24_rep3",
  "TG_p24_rep1",
  "TG_p24_rep2",
  "TG_p24_rep3",
  "Control_p29_rep1",
  "Control_p29_rep2",
  "Control_p29_rep3",
  "TG_p29_rep1",
  "TG_p29_rep2",
  "TG_p29_rep3"
)

# sample_metadata.dt$alias%in%colnames(counts.mtx)
# counts.mtx <- counts.mtx[,sample_metadata.dt$sample]

#################################
## Create SingleCellExperiment ##
#################################

coldata_to_sce <- DataFrame(sample_metadata.dt[alias%in%colnames(counts.mtx)] %>% setkey(alias) %>% .[colnames(counts.mtx)] %>% tibble::column_to_rownames("alias"))
stopifnot(colnames(counts.mtx)==rownames(coldata_to_sce))

rowdata_to_sce <- DataFrame(gene_metadata.dt[transcript_id%in%rownames(counts.mtx)] %>% setkey(transcript_id) %>% .[rownames(counts.mtx)]  %>% tibble::column_to_rownames("transcript_id"))
stopifnot(rownames(counts.mtx)==rownames(rowdata_to_sce))

sce <- SingleCellExperiment(
  assays = list(tpm=counts.mtx), 
  rowData = rowdata_to_sce,
  colData = coldata_to_sce
)

##########
## Save ##
##########

saveRDS(sce, io$outfile)
