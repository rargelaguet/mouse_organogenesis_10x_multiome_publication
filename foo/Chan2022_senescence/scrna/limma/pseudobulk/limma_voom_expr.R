source(here::here("settings.R"))
source(here::here("utils.R"))

library(limma)

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"results/cell_cycle/sample_metadata_with_cell_cycle_scores.txt")
io$sce <- file.path(io$basedir,"rna/processed/salmon/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/limma_voom"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$cell_lines <- c("HC1", "HC2", "HC3")
opts$treatments <- c("Control")
opts$percent_oxygen <- c(21)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata) %>%
  .[treatment%in%opts$treatments & cell_line%in%opts$cell_lines & percent_oxygen%in%opts$percent_oxygen]

#########################
## Load gene expression ##
#########################

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]
# rowData(sce)$ens_id <- rownames(sce)
# rownames(sce) <- rowData(sce)$gene_name
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

######################################
## limma voom, regressing out donor ##
######################################

# Create a design matrix
design <- model.matrix(~ passage, data = colData(sce))

# Feature selection
keep <- edgeR::filterByExpr(sce, design)
sce.filt <- sce[keep,]

# voom
v <- voom(assay(sce.filt,"counts"), design = design, plot = TRUE)

# Add log cpm to the SingleCellExperiment
assay(sce.filt,"log_cpm") <- v$E

# Estimate the inter-media correlation
corfit <- duplicateCorrelation(v$E, design = design, block = colData(sce)$cell_line, weights = v$weights)

# Fit the linear model with the estimated correlation and voom weights
fit_mixed <- lmFit(v$E, design, block = colData(sce)$cell_line, correlation = corfit$consensus.correlation, weights = v$weights)

# Perform empirical Bayes moderation
fit_mixed <- eBayes(fit_mixed)

# Extract the coefficients and p-values
# (The rest of your result extraction code remains the same)
limma_results_list <- list()
limma_results_list[["passage"]] <- topTable(fit_mixed, coef='passage', number=Inf)[,c("logFC","P.Value","adj.P.Val")] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("feature","coef","pval","pval_adj"))

#####################################################################
## limma voom, including proliferation score, regressing out donor ##
#####################################################################

# Create a design matrix
design <- model.matrix(~ passage + cell_cycle_score, data = colData(sce))

# Feature selection
keep <- edgeR::filterByExpr(sce, design)
sce.filt <- sce[keep,]

# voom
v <- voom(assay(sce.filt,"counts"), design = design, plot = TRUE)

# Add log cpm to the SingleCellExperiment
assay(sce.filt,"log_cpm") <- v$E

# Estimate the inter-media correlation
corfit <- duplicateCorrelation(v$E, design = design, block = colData(sce)$cell_line, weights = v$weights)

# Fit the linear model with the estimated correlation and voom weights
fit_mixed <- lmFit(v$E, design, block = colData(sce)$cell_line, correlation = corfit$consensus.correlation, weights = v$weights)

# Perform empirical Bayes moderation
fit_mixed <- eBayes(fit_mixed)

# Extract the coefficients and p-values
# (The rest of your result extraction code remains the same)
limma_results_proliferation_list <- list()
limma_results_proliferation_list[["passage"]] <- topTable(fit_mixed, coef='passage', number=Inf)[,c("logFC","P.Value","adj.P.Val")] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("feature","coef","pval","pval_adj"))
limma_results_proliferation_list[["cell_cycle"]] <- topTable(fit_mixed, coef='cell_cycle_score', number=Inf)[,c("logFC","P.Value","adj.P.Val")] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("feature","coef","pval","pval_adj"))

#################################
## update SingleCellExperiment ##
#################################

metadata(sce.filt)[["C1__passage"]] <- limma_results_list[["passage"]]
metadata(sce.filt)[["C1__proliferation-corrected-passage"]] <- limma_results_proliferation_list[["passage"]]

metadata(sce.filt)[["analysis_descriptions"]] <- list(
  "C1__passage" = "C1__passage: Limma voom model fit to passage number, regressing out donor variability",
  "C1__proliferation-corrected-passage" = "Limma voom model fit to passage number, regressing out donor variability and including proliferation score as a covariate."
)

##########
## Save ##
##########
  
fwrite(limma_results_list[["passage"]], file.path(io$outdir,"limma_results_passage.txt.gz"), sep="\t", quote=F, na="NA")
fwrite(limma_results_proliferation_list[["passage"]], file.path(io$outdir,"limma_results_passage_proliferation_corrected.txt.gz"), sep="\t", quote=F, na="NA")
saveRDS(sce.filt, file.path(io$outdir,"SingleCellExperiment_limma_voom.rds"))

