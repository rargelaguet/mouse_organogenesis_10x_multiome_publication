library(scran)
library(variancePartition)

source(here::here("settings.R"))
source(here::here("utils.R"))

set.seed(42)

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"processed/rna/sample_metadata.txt")
io$sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/rna/variance_partition"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Parallel processing
# cl <- makeCluster(12)
# registerDoParallel(cl)
BPPARAM <- BiocParallel::MulticoreParam(workers=8)

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata)

#########################
## Load RNA expression ##
#########################

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]
# colnames(sce) <- sample_metadata.dt$sample
# rownames(sce) <- rowData(sce)$symbol

colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

###################
## Normalisation ##
###################

assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))+1) # log cpm

#######################
## Feature selection ##
#######################

opts$ngenes <- 100

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=opts$ngenes) %>% rownames
sce_filt <- sce[hvgs,]

########################
## Variance partition ##
########################

# Define equation
eq <- ~ (1|Treatment) + Passage

# Define metadata
info <- sample_metadata.dt[,c("Treatment","Passage")] %>% as.data.frame
rownames(info) <- sample_metadata.dt$sample

# Run
varPart <- fitExtractVarPartModel(assay(sce_filt,"logcounts"), formula=eq, data=info, BPPARAM=BPPARAM)

##########
## Save ##
##########

saveRDS(varPart, file.path(io$outdir,"variance_partition.rds"))

