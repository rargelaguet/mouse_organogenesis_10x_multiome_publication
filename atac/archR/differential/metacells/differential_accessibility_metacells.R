here::i_am("atac/archR/differential/metacells/differential_accessibility_metacells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(edgeR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atac_matrix_file',        type="character",                               help='')
p$add_argument('--metadata',        type="character",                               help='')
p$add_argument('--matrix',          type="character",  default="PeakMatrix",   help='Matrix to use')
p$add_argument('--group_variable',          type="character",   help='Group variable')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--celltypes',    type="character",    default="all", nargs="+", help='Celltypes to use')
p$add_argument('--samples',    type="character",    default="all", nargs="+", help='Samples to use')
p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per group')
p$add_argument('--outfile',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$matrix <- "PeakMatrix" # "GeneScoreMatrix_TSS"
# args$atac_matrix_file <- file.path(io$basedir,sprintf("results/atac/archR/metacells/all_cells/PeakMatrix/%s_summarized_experiment_metacells.rds",args$matrix))
# args$metadata <- file.path(io$basedir,sprintf("results/atac/archR/metacells/all_cells/PeakMatrix/metacells_metadata.txt.gz"))
# # args$group_variable <- "genotype"; args$groupA <- "WT"; args$groupB <- "T_KO"
# args$group_variable <- "celltype"; args$groupA <- "Neural_crest"; args$groupB <- "Epiblast"
# args$samples <- "all" # c("E8.5_CRISPR_T_KO", "E8.5_CRISPR_T_WT")
# args$celltypes <- "all"
# args$min_cells <- 5
# args$outfile <- file.path(io$basedir,sprintf("results/atac/archR/differential/metacells/%s/%s/%s_vs_%s.txt.gz",args$group_variable,args$matrix,args$groupA,args$groupB))
## END TEST ##

#####################
## Parse arguments ##
#####################

# Define cell types
if (args$celltypes[1]=="all") {
  args$celltypes <- opts$celltypes
} else {
  # stopifnot(args$celltypes%in%c(opts$celltypes,"Erythroid","Blood_progenitors"))
  stopifnot(args$celltypes%in%opts$celltypes)
}

# Define samples
if (args$samples[1]=="all") {
  args$samples <- opts$samples
} else {
  stopifnot(args$samples%in%opts$samples)
}

dir.create(dirname(args$outfile), showWarnings=F, recursive = T)

# stupid stuff but otherwise the snakemake pipeline doesn't work
if (args$groupA==args$groupB) {
  warning("groupA and groupB are the same, saving an empty file...")
  out <- data.table(feature=NA, logFC=NA, padj_fdr=NA, groupA_N=NA, groupB_N=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  quit(status=0)
}

#####################
## Define settings ##
#####################

# Options
opts$groups <- c(args$groupA,args$groupB)
opts$min_cdr <- 0.30

# Load utils
source(here::here("atac/archR/differential/metacells/utils.R"))

###################
## Load metadata ##
###################

print("Loading metacells...")

metacell_metadata.dt <- fread(args$metadata) %>%
  .[celltype%in%args$celltypes & sample%in%args$samples] %>%
  .[,celltype_genotype:=sprintf("%s_%s",celltype,genotype)]

stopifnot(args$group_variable%in%colnames(metacell_metadata.dt))

metacell_metadata.dt <- metacell_metadata.dt %>%
  .[,group:=eval(as.name(args$group_variable))] %>%
  .[group%in%opts$groups] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

# print(table(metacell_metadata.dt$sample))
# print(table(metacell_metadata.dt$celltype))
# print(table(metacell_metadata.dt$group))

tmp <- table(metacell_metadata.dt$group)
if (any(tmp<=args$min_cells)) {
  warning("Not enough cells per group to perform DE, saving an empty file...")
  out <- data.table(feature=NA, logFC=NA, padj_fdr=NA, groupA_N=tmp[args$groupA], groupB_N=tmp[args$groupB])
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  quit(status=0)
} 

######################
## Load ATAC Matrix ##
######################

print(sprintf("Fetching ATAC matrix: '%s'...",args$atac_matrix_file))

# Load 
atac.se <- readRDS(args$atac_matrix_file)[,metacell_metadata.dt$metacell]
colData(atac.se) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

# Normalise
assayNames(atac.se)[1] <- "counts"
assay(atac.se,"logcounts") <- log(1e6*(sweep(assay(atac.se,"counts"),2,colSums(assay(atac.se,"counts")),"/"))+1)

##########################################
## calculate detection rate per feature ##
##########################################

acc.dt <- data.table(
  feature = rownames(atac.se),
  detection_rate_groupA = rowMeans(assay(atac.se,"counts")[,atac.se$celltype==args$groupA]>0) %>% round(2),
  detection_rate_groupB = rowMeans(assay(atac.se,"counts")[,atac.se$celltype==args$groupB]>0) %>% round(2),
  mean_groupA = rowMeans(assay(atac.se[,atac.se$celltype==args$groupA],"logcounts")) %>% round(2),
  mean_groupB = rowMeans(assay(atac.se[,atac.se$celltype==args$groupB],"logcounts")) %>% round(2)
)

###################################################
## Differential accessibility testing with edgeR ##
##################################################

# Consider only features that have a minimum detection rate
features.to.use <- acc.dt[detection_rate_groupA>=opts$min_cdr | detection_rate_groupB>=opts$min_cdr,feature]

# Convert SCE to DGEList
atac_metacells.dge <- DGEList(assay(atac.se[features.to.use,],"counts"))
# sce_edger <- scran::convertTo(atac.se, type="edgeR")

# Define design matrix (with intercept)
cdr <- colMeans(assay(atac.se,"counts")>0)
design <- model.matrix(~cdr+atac.se$group)

# Estimate dispersions
atac_metacells.dge  <- estimateDisp(atac_metacells.dge,design)

# Fit GLM
fit <- glmQLFit(atac_metacells.dge,design)

# Likelihood ratio test
lrt <- glmQLFTest(fit)

# Construct output data.frame
out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
  setnames(c("feature","logFC","logCPM","LR","p.value","padj_fdr")) %>%
  .[,c("logCPM","LR","p.value"):=NULL] %>%
  .[,c("groupA_N","groupB_N"):=list(table(metacell_metadata.dt$group)[1],table(metacell_metadata.dt$group)[2])]%>% 
  .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3),round(logFC,3))] %>%
  merge(acc.dt[,c("feature","mean_groupA","mean_groupB")], by="feature", all.y=TRUE) %>%
  .[is.na(logFC),c("logFC","padj_fdr"):=list(0,1)] %>%
  setorder(padj_fdr, na.last=T)


##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)


##########
## TEST ##
##########

# i <- "chr6:91880813-91881413"

# out[feature==i]
# tmp <- atac.se[]
# mean(log(assay(tmp[,tmp$genotype=="WT"])+1))
# mean(log(assay(tmp[,tmp$genotype=="T_KO"])+1))

# atac_metacells.dge$AveLogCPM
