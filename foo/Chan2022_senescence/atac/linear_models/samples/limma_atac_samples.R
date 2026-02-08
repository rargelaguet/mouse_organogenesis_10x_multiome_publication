source(here::here("settings.R"))
source(here::here("utils.R"))

library(limma)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
# io$sce <- paste0(io$basedir,"/processed/atac/SingleCellExperiment.rds")
io$sce <- paste0(io$basedir,"/processed/atac/SummarizedExperiment_ATAC_Jack.rds")
# io$outdir <- file.path(io$basedir,"results/atac/limma/all_samples"); dir.create(io$outdir, showWarnings = F, recursive = T)
io$outdir <- file.path(io$basedir,"results/atac/limma/samples/control_samples"); dir.create(io$outdir, showWarnings = F, recursive = T)

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata) %>%
  .[Treatment=="Control"]

###############
## Load ATAC ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]
# colnames(sce) <- sample_metadata.dt$sample
tmp <- as.data.table(rowRanges(sce)) %>% .[,c("seqnames","start","end","name")] %>% .[,feature:=sprintf("%s:%s-%s",seqnames,start,end)]
rownames(sce) <- tmp$feature

colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

assay(sce,"counts") <- assay(sce,"rawCounts")
###########################
## Load feature metadata ##
###########################

# gene_metadata.dt <- fread(io$gene_metadata)

# gene_metadata_champ.dt <- fread("/Users/rargelaguet/data/Chan2022_senescence/processed/rna/feature_metadata.txt.gz")

###################
## Normalisation ##
###################

assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(assay(sce,"counts")),"/"))+1) # log cpm
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

atac.mtx <- assay(sce,"logcounts")
# hist(atac.mtx)

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

# View(limma_results[["Passage:Treatment"]])
# limma_results[["Passage"]][feature=="SOX11"]

##########
## Save ##
##########

saveRDS(limma_results, file.path(io$outdir,"limma_results.rds"))

