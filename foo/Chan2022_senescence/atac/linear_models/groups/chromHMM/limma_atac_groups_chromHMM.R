source(here::here("settings.R"))
source(here::here("utils.R"))

library(limma)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
io$feature_metadata <- "/Users/rargelaguet/data/hg38_regulation/chromHMM/HDF_6_dense.merged.bed"
io$sce <- file.path(io$basedir,"results/atac/deeptools/SingleCellExperiment_chromHMM.rds")
io$outdir <- file.path(io$basedir,"results/atac/limma/groups/control_groups"); dir.create(io$outdir, showWarnings = F, recursive = T)

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata) %>%
  .[Group%in%opts$groups,c("Group","Passage","Treatment")] %>% 
  .[,group:=sprintf("%s_p%s",Treatment,Passage)] %>%
  unique

sample_metadata.dt <- sample_metadata.dt[Treatment=="Control"]

###############
## Load ATAC ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$group]

colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("group") %>% DataFrame

###########################
## Load feature metadata ##
###########################

feature_metadata.dt <- fread(io$feature_metadata, select=c(1,2,3,4)) %>%
  setnames(c("chr","start","end","annotation")) %>%
  .[,feature:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("feature","annotation","chr","start","end")] %>%
  .[!feature%in%duplicated.sites] %>% setkey(feature) %>% .[rownames(sce)]

###################
## Normalisation ##
###################

# assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))+1) # log cpm
# assay(sce,"logcounts") <- 1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/")) # cpm

#################################################
## Fit linear model with no interaction effect ##
#################################################

# # genes.to.plot <- which(apply(assay(sce,"logcounts"),1,mean)>=6) %>% names
# genes.to.plot <- rownames(sce)
# 
# # dt <- fread(file.path(io$outdir,"linear_models_te_age.tsv.gz"))
# 
# donors.to.plot <- unique(sample_metadata.dt$group)
# 
# # i <- "RMER13A1:ERVK:LTR"; j <- "Brain"
# # i <- "L1"; j <- "Brain"
# dt <- genes.to.plot %>% map(function(i) {
#   donors.to.plot %>% map(function(j) {
#     print(sprintf("%s - %s",i,j))
#     
#     to.plot <- data.table(
#       sample = colnames(sce),
#       expr = assay(sce,"logcounts")[i,]
#     ) %>% merge(sample_metadata.dt[group==j,c("sample","passage","group")],by="sample") 
#     
#     model <- lm(expr~passage, data=to.plot)
#   
#     data.table(
#       feature = i,
#       group = j,
#       coef = summary(model)$coefficients["passage","Estimate"],
#       p = summary(model)$coefficients["passage","Pr(>|t|)"],
#       r = cor(to.plot$expr, to.plot$passage)
#     )
#   }) %>% rbindlist
# }) %>% rbindlist

###########
## limma ##
###########

# atac.mtx <- assay(sce,"logcounts")
atac.mtx <- assay(sce)

# Create a design matrix
# design <- model.matrix(~ Passage + Treatment + Passage:Treatment, data = colData(sce))
design <- model.matrix(~ Passage, data = colData(sce))

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

limma_results_passage.dt <- limma_results[["Passage"]] %>% as.data.table %>%
  merge(feature_metadata.dt, by="feature")

##########
## Save ##
##########

saveRDS(limma_results, file.path(io$outdir,"limma_results.rds"))

