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
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# io$basedir <- file.path(io$basedir,"test")
# args$matrix <- "PeakMatrix"
# args$atac_matrix_file <- file.path(io$basedir,sprintf("results/atac/archR/pseudobulk/celltype/%s/%s_pseudobulk_with_replicates.rds",args$matrix,args$matrix))
# args$groupA <- "Epiblast"
# args$groupB <- "Neural_crest"
# args$outfile <- NULL
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
  out <- data.table(feature=NA, logFC=NA, padj_fdr=NA, mean_groupA=NA, mean_groupB=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  warning("groupA and groupB are the same, saving an empty file...")
  quit(status=0)
}

###############
## Load data ##
###############

print(sprintf("Fetching ATAC matrix: '%s'...",args$atac_matrix_file))

atac.se <- readRDS(args$atac_matrix_file)

# parse
if (!"celltype" %in% colnames(colData(atac.se))) {
  atac.se$celltype <- colnames(atac.se) %>% strsplit("-") %>% map_chr(1)
}
atac.se <- atac.se[,atac.se$celltype%in%opts$groups]
atac.se$celltype <- factor(atac.se$celltype, levels=opts$groups)

# Normalise with log2 counts (for consistentcy with edgeR)
assayNames(atac.se)[1] <- "counts"
logcounts(atac.se) <- log2(1e6*(sweep(counts(atac.se),2,colSums(counts(atac.se)),"/"))+1)

# check that we have pseudobulk replicates for both cell types
if (any(!opts$groups%in%unique(atac.se$celltype))) {
  out <- data.table(feature=NA, logFC=NA, padj_fdr=NA, mean_groupA=NA, mean_groupB=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  warning("groups not found, saving an empty file...")
  quit(status=0)
}

##################################
## calculate feature statistics ##
##################################

acc.dt <- data.table(
  feature = rownames(atac.se),
  detection_rate_groupA = rowMeans(assay(atac.se,"counts")[,atac.se$celltype==opts$groups[1]]>0) %>% round(2),
  detection_rate_groupB = rowMeans(assay(atac.se,"counts")[,atac.se$celltype==opts$groups[2]]>0) %>% round(2),
  mean_groupA = rowMeans(logcounts(atac.se[,atac.se$celltype==args$groupA])) %>% round(2),
  mean_groupB = rowMeans(logcounts(atac.se[,atac.se$celltype==args$groupB])) %>% round(2)
)

#####################################
## Differential testing with edgeR ##
#####################################

# Consider only features that have a minimum detection rate
features.to.use <- acc.dt[detection_rate_groupA>=opts$min_cdr | detection_rate_groupB>=opts$min_cdr,feature]

# Convert SCE to DGEList
atac_metacells.dge <- DGEList(assay(atac.se[features.to.use,],"counts"))

# Define design matrix (with intercept)
design <- model.matrix(~atac.se$celltype)

# Estimate dispersions
atac_metacells.dge <- estimateDisp(atac_metacells.dge,design)

# Fit GLM
fit <- glmQLFit(atac_metacells.dge,design)

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
