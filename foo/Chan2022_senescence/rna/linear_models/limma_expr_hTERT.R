source(here::here("settings.R"))
source(here::here("utils.R"))

library(limma)

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"rna/results/limma"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$plot_logcounts <- TRUE

opts$conditions <- c(
  "hTERT"
)

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] 

sample_metadata.dt[,time_numeric:=as.numeric(gsub("TP","",time))] 

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

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

if (opts$plot_logcounts) {
  assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))+1) # log cpm
} else { 
  assay(sce,"logcounts") <- 1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/")) # cpm
}

#######################
## Feature selection ##
#######################

tmp <- apply(assay(sce,"logcounts"),1,var,na.rm=T)
sce.filt <- sce[tmp>0.01,]

###########
## limma ##
###########

expr.mtx <- assay(sce.filt,"logcounts")

# Create a design matrix
design <- model.matrix(~ time_numeric, data = colData(sce.filt))

# Fit the linear model
fit <- lmFit(expr.mtx, design)

# Perform empirical Bayes moderation
fit <- eBayes(fit)

# Extract the coefficients and p-values
limma_results.dt <- topTable(fit, coef='time_numeric', number=Inf)[,c("logFC","P.Value","adj.P.Val")] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("ens_id","coef","pval","pval_adj")) %>%
  merge(gene_metadata.dt[,c("ens_id","symbol")]) %>% 
  setorder(pval)

##########
## Save ##
##########

fwrite(limma_results.dt, file.path(io$outdir,"limma_results_hTERT.txt.gz"), sep="\t", quote=F, na="NA")

