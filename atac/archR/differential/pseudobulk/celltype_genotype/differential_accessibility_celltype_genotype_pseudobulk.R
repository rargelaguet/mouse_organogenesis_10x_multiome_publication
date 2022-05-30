here::i_am("atac/archR/differential/pseudobulk/celltype/differential_accessibility_pseudobulk.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(edgeR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--matrix',        type="character",                               help='')
p$add_argument('--atac_matrix_file',        type="character",                               help='')
p$add_argument('--celltypes',    type="character",  nargs="+",  help='celltype')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# io$basedir <- file.path(io$basedir,"test")
# args$matrix <- "PeakMatrix"
# args$atac_matrix_file <- file.path(io$basedir,sprintf("results/atac/archR/pseudobulk/celltype_genotype/%s/%s_pseudobulk_with_replicates.rds",args$matrix,args$matrix))
# args$celltypes <- "Neural_crest"; args$groupA <- "WT"; args$groupB <- "T_KO"
# args$outfile <- tempfile()
## END TEST

dir.create(dirname(args$outfile), showWarnings=F, recursive = T)

#####################
## Define settings ##
#####################

# Define groups
opts$groups <- c(args$groupA,args$groupB)

opts$min_cdr <- 0.30

# stupid stuff but otherwise the snakemake pipeline doesn't work
if (args$groupA==args$groupB) {
  out <- data.table(feature=NA, logFC=NA, padj_fdr=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  warning("groupA and groupB are the same, saving an empty file...")
  quit(status=0)
}

###############
## Load data ##
###############

print(sprintf("Fetching ATAC matrix: '%s'...",args$atac_matrix_file))

atac.se <- readRDS(args$atac_matrix_file)
assayNames(atac.se)[1] <- "counts"

# Normalise with log2 counts (for consistentcy with edgeR)
logcounts(atac.se) <- log2(1e6*(sweep(counts(atac.se),2,colSums(counts(atac.se)),"/"))+1)

# parse
if (!"celltype" %in% colnames(colData(atac.se))) {
  atac.se$celltype <- colnames(atac.se) %>% strsplit("-") %>% map_chr(1)
}
if (!"genotype" %in% colnames(colData(atac.se))) {
  atac.se$genotype <- colnames(atac.se) %>% strsplit("-") %>% map_chr(2) %>% strsplit("_rep") %>% map_chr(1)
}
if (!"celltype_genotype" %in% colnames(colData(atac.se))) {
  atac.se$celltype_genotype <- sprintf("%s_%s",atac.se$celltype,atac.se$genotype)
}

# temporary (otherwise subsetting doesn't work???)
if ("start" %in% colnames(rowData(atac.se))) {
  rowData(atac.se)$start <- NULL
}
if ("end" %in% colnames(rowData(atac.se))) {
  rowData(atac.se)$end <- NULL
}

###################
## Sanity checks ##
###################

tmp <- sprintf("%s_%s",args$celltypes,opts$groups)
if (any(!tmp%in%unique(atac.se$celltype_genotype))) {
  warning("...")
  out <- data.table(feature=NA, logFC=NA, padj_fdr=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  quit(status=0)
} else {
  atac.se <- atac.se[,atac.se$celltype %in% args$celltypes]
}

# # subset genotypes
stopifnot(opts$groups%in%atac.se$genotype)
atac.se$genotype <- factor(atac.se$genotype, levels=opts$groups)
print(table(atac.se$celltype_genotype))

##################################
## calculate feature statistics ##
##################################

acc.dt <- data.table(
  feature = rownames(atac.se),
  detection_rate_groupA = rowMeans(assay(atac.se,"counts")[,atac.se$genotype==opts$groups[1]]>0) %>% round(2),
  detection_rate_groupB = rowMeans(assay(atac.se,"counts")[,atac.se$genotype==opts$groups[2]]>0) %>% round(2),
  mean_groupA = rowMeans(logcounts(atac.se[,atac.se$genotype==args$groupA])) %>% round(2),
  mean_groupB = rowMeans(logcounts(atac.se[,atac.se$genotype==args$groupB])) %>% round(2)
)

#####################################
## Differential testing with edgeR ##
#####################################

# Consider only features that have a minimum detection rate
features.to.use <- acc.dt[detection_rate_groupA>=opts$min_cdr | detection_rate_groupB>=opts$min_cdr,feature]

# Create DGEList
atac.dge <- DGEList(assay(atac.se[features.to.use,],"counts"))

# Define design matrix (with intercept)
design <- model.matrix(~atac.se$genotype)

# Estimate dispersions
atac.dge <- estimateDisp(atac.dge,design)

# Fit GLM
fit <- glmQLFit(atac.dge,design)

# Likelihood ratio test
lrt <- glmQLFTest(fit)

# Construct output data.frame
out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
	setnames(c("feature","logFC","logCPM","LR","p.value","padj_fdr")) %>%
	.[,c("logCPM","LR","p.value"):=NULL] %>%
  merge(acc.dt[,c("feature","mean_groupA","mean_groupB")], by="feature", all.y=TRUE) %>%
	.[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3),round(logFC,3))] %>%
	.[is.na(logFC),c("logFC","padj_fdr"):=list(0,1)] %>%
	setorder(padj_fdr, na.last=T)

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
