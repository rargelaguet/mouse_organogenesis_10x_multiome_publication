here::i_am("rna/differential/pseudobulk/celltype_genotype/differential_celltype_genotype_pseudobulk.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(edgeR))
suppressMessages(library(scater))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',    type="character",    help='SingleCellExperiment file')
p$add_argument('--celltypes',    type="character",  nargs="+",  help='celltype')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# io$basedir <- file.path(io$basedir,"test")
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype_genotype/SingleCellExperiment_pseudobulk_with_replicates.rds")
# args$celltypes <- c("Intermediate_mesoderm")
# args$groupA <- "WT"
# args$groupB <- "T_KO"
# args$outfile <- tempfile()
## END TEST

dir.create(dirname(args$outfile), showWarnings = F)

# Define groups
opts$groups <- c(args$groupA,args$groupB)

#########################
## Load RNA expression ##
#########################

# Load SingleCellExperiment object
sce <- readRDS(args$sce)

# parse annotations
if (!"celltype" %in% colnames(colData(sce))) {
  sce$celltype <- colnames(sce) %>% strsplit("-") %>% map_chr(1)
}
if (!"genotype" %in% colnames(colData(sce))) {
  sce$genotype <- colnames(sce) %>% strsplit("-") %>% map_chr(2) %>% strsplit("_rep") %>% map_chr(1)
}
if (!"celltype_genotype" %in% colnames(colData(sce))) {
  sce$celltype_genotype <- sprintf("%s_%s",sce$celltype,sce$genotype)
}

###################
## Sanity checks ##
###################

tmp <- sprintf("%s_%s",args$celltypes,opts$groups)
if (any(!tmp%in%unique(sce$celltype_genotype))) {
  warning("...")
  out <- data.table(feature=NA, logFC=NA, padj_fdr=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  quit(status=0)
} else {
  sce <- sce[,sce$celltype %in% args$celltypes]
}

# subset genotypes
stopifnot(opts$groups%in%sce$genotype)
sce$genotype <- factor(sce$genotype, levels=opts$groups)
print(table(sce$celltype_genotype))

#########################################
## Calculate average expression levels ##
#########################################

expr.dt <- data.table(
  gene = rownames(sce),
  mean_groupA = rowMeans(logcounts(sce[,sce$genotype==args$groupA])) %>% round(2),
  mean_groupB = rowMeans(logcounts(sce[,sce$genotype==args$groupB])) %>% round(2)
)

#######################
## Feature selection ##
#######################

opts$min.expr <- 4 # 2**4 = 16, at least an average of 8 counts per milion for each group

genes.to.use <- expr.dt[mean_groupA>=opts$min.expr | mean_groupB>=opts$min.expr,gene]

################################################
## Differential expression testing with edgeR ##
################################################

# Convert SCE to DGEList
sce_edger <- scran::convertTo(sce[genes.to.use,], type="edgeR")

# Define design matrix (with intercept)
design <- model.matrix(~sce$genotype)

# Estimate dispersions
sce_edger  <- estimateDisp(sce_edger,design)

# Fit GLM
fit <- glmQLFit(sce_edger,design)

# Likelihood ratio test
lrt <- glmQLFTest(fit)

# Construct output data.frame
out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
  setnames(c("gene","logFC","logCPM","LR","p.value","padj_fdr")) %>%
  .[,c("logCPM","LR","p.value"):=NULL] %>%
  .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3), round(logFC,3))] %>%
  merge(expr.dt, by="gene", all.y=TRUE) %>%
  setorder(padj_fdr, na.last=T)

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
