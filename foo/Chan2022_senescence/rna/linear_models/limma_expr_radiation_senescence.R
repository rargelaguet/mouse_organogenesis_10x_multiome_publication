source(here::here("settings.R"))
source(here::here("utils.R"))

library(limma)

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"rna/results/limma"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$conditions <- c(
  "Radiation_Senescence"
)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] 

sample_metadata.dt[time=="D2_ctrl", time:="D0"]
sample_metadata.dt[,time_numeric:=as.numeric(gsub("D","",time))] 

table(sample_metadata.dt$time_numeric)

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

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata) %>%
  .[ens_id%in%rownames(sce)]

###################
## Normalisation ##
###################

assay(sce,"cpm") <- 1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))
assay(sce,"log_cpm") <- log2(assay(sce,"cpm")+1)

#######################
## Feature selection ##
#######################

tmp <- apply(assay(sce,"log_cpm"),1,var,na.rm=T)
sce.filt <- sce[tmp>0.01,]

###########
## limma ##
###########

expr.mtx <- assay(sce.filt,"log_cpm")

# Create a design matrix
design <- model.matrix(~ time_numeric, data = colData(sce.filt))

# Fit the linear model
fit <- lmFit(expr.mtx, design)

# Perform empirical Bayes moderation
fit <- eBayes(fit)

# Extract the coefficients and p-values
limma_results.dt <- topTable(fit, coef='time_numeric', number=Inf)[,c("logFC","P.Value","adj.P.Val")] %>% as.data.table(keep.rownames = T) %>%
  setnames(c("ens_id","coef","pval","pval_adj")) %>%
  merge(gene_metadata.dt[,c("ens_id","symbol")], all.y=TRUE) %>% 
  setorder(pval)

##########
## Save ##
##########

fwrite(limma_results.dt, file.path(io$outdir,"limma_results_radiation_senescence.txt.gz"), sep="\t", quote=F, na="NA")

#####################
## Save (extended) ##
#####################

limma_results_extended.dt <- topTable(fit, coef = "time_numeric", number = Inf) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn", "ens_id")

gene_metadata_tmp.dt <- gene_metadata.dt[, c("chr", "ens_id", "description", "symbol")] %>%
  setnames("symbol", "gene_name")

limma_results_extended_save.dt <- merge(limma_results_extended.dt, gene_metadata_tmp.dt, by = "ens_id") %>%
  .[, c("gene_name", "ens_id", "chr", "description", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

fwrite(limma_results_extended_save.dt, file.path(io$outdir,"extended/limma_results_radiation_senescence_extended.txt.gz"), sep="\t", quote=F, na="NA")
