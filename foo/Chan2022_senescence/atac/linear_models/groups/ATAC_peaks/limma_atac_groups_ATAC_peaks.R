source(here::here("settings.R"))
source(here::here("utils.R"))

library(limma)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
io$feature_metadata <- file.path(io$basedir,"processed/atac/jack/feature_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/atac/jack/SingleCellExperiment.rds")
io$outfile <- file.path(io$basedir,"results/atac/limma/samples/all_samples/limma_results_ATAC_peaks.rds")

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata) %>%
  .[Group%in%opts$groups] %>% 
  .[,group:=sprintf("%s_p%s",Treatment,Passage)] %>%
  unique

# sample_metadata.dt <- sample_metadata.dt[Treatment=="Control"]

###############
## Load ATAC ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]

colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

###########################
## Load feature metadata ##
###########################

# feature_metadata.dt <- fread(io$feature_metadata) %>%
#   setnames(c("chr","start","end","annotation")) %>%
#   .[,feature:=sprintf("%s:%s-%s",chr,start,end)] %>%
#   .[,c("feature","annotation","chr","start","end")] %>%
#   .[!feature%in%duplicated.sites] %>% setkey(feature) %>% .[rownames(sce)]

###################
## Normalisation ##
###################

assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))+1) # log cpm
# assay(sce,"logcounts") <- 1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/")) # cpm

###########
## limma ##
###########

atac.mtx <- assay(sce,"logcounts")
# atac.mtx <- assay(sce)

# Create a design matrix
design <- model.matrix(~ Passage + Treatment + Passage:Treatment, data = colData(sce))
# design <- model.matrix(~ Passage, data = colData(sce))

# Fit the linear model
fit <- lmFit(atac.mtx, design)

# Perform empirical Bayes moderation
fit <- eBayes(fit)

# Extract the coefficients and p-values
# head(fit$coefficients)
limma_results <- list()
limma_results[["Passage"]] <- topTable(fit, coef='Passage', number=Inf)[,c("logFC","adj.P.Val")] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("feature","logFC","pval_adj"))
limma_results[["Treatment"]] <- topTable(fit, coef='TreatmentTG', number=Inf)[,c("logFC","adj.P.Val")] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("feature","logFC","pval_adj"))
limma_results[["Passage:Treatment"]] <- topTable(fit, coef='Passage:TreatmentTG', number=Inf)[,c("logFC","adj.P.Val")] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("feature","logFC","pval_adj"))

#######################
## Create data.table ##
#######################

limma_results_passage.dt <- limma_results[["Passage"]]# %>% merge(feature_metadata.dt, by="feature")

##########
## Save ##
##########

# saveRDS(limma_results, file.path(io$outdir,"limma_results_BJ_reprogramming.rds"))
saveRDS(limma_results, io$outfile)

