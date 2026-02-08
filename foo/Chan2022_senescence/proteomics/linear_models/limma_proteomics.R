source(here::here("settings.R"))
source(here::here("utils.R"))

library(limma)

#####################
## Define settings ##
#####################

# I/O
io$sample_metadata <- file.path(io$basedir,"processed/proteomics/sample_metadata.txt")
io$feature_metadata <- file.path(io$basedir,"processed/proteomics/feature_metadata.txt")
io$sce <- paste0(io$basedir,"/processed/proteomics/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/proteomics/limma"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$plot_logcounts <- TRUE

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$sample_metadata)

###############
## Load data ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$alias]
# colnames(sce) <- sample_metadata.dt$sample
rownames(sce) <- rowData(sce)$protein_entry

# colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

###########################
## Load feature metadata ##
###########################

feature_metadata.dt <- fread(io$feature_metadata)

################
## Imputation ##
################

# Remove features with large fraction of NAs
tmp <- rowMeans(is.na(assay(sce,"proteins_normalised")))
sce <- sce[tmp<0.25]

# imputation
# assay(sce,"proteins_normalised_imputed") <- impute::impute.knn(assay(sce,"proteins_normalised"), k=10, rng.seed=42)
assay(sce,"proteins_normalised_imputed") <- impute::impute.knn(assay(sce,"proteins_normalised"), k=3, rng.seed=42)[[1]]

###################
## Normalisation ##
###################

# library(matrixStats)

# z-scores
assay(sce,"zscores") <- (assay(sce,"proteins_normalised")-rowMeans(assay(sce,"proteins_normalised"),na.rm=T)) / sqrt(rowVars(assay(sce,"proteins_normalised"),na.rm=T))


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

expr.mtx <- assay(sce,"proteins_normalised")

# Create a design matrix
design <- model.matrix(~ Passage + Treatment + Passage:Treatment, data = colData(sce))

# Fit the linear model
fit <- lmFit(expr.mtx, design)

# Perform empirical Bayes moderation
fit <- eBayes(fit)

# Extract the coefficients and p-values
# head(fit$coefficients)
limma_results <- list()
limma_results[["Passage"]] <- topTable(fit, coef='Passage', number=Inf)[,c("ID","logFC","adj.P.Val")] %>% as.data.table(keep.rownames = F) %>%
  setnames(c("feature","logFC","pval_adj"))
limma_results[["Treatment"]] <- topTable(fit, coef='TreatmentTG', number=Inf)[,c("ID","logFC","adj.P.Val")] %>% as.data.table(keep.rownames = F) %>%
  setnames(c("feature","logFC","pval_adj"))
limma_results[["Passage:Treatment"]] <- topTable(fit, coef='Passage:TreatmentTG', number=Inf)[,c("ID","logFC","adj.P.Val")] %>% as.data.table(keep.rownames = F) %>%
  setnames(c("feature","logFC","pval_adj"))

# View(limma_results[["Passage:Treatment"]])
# limma_results[["Passage"]][feature=="SOX11"]

##########
## Save ##
##########

saveRDS(limma_results, file.path(io$outdir,"limma_results.rds"))

