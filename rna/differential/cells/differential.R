here::i_am("rna/differential/cells/differential.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(edgeR))
suppressMessages(library(scater))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='Cell metadata file')
p$add_argument('--sce',    type="character",    help='SingleCellExperiment file')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--celltypes',    type="character",    default="all", nargs="+", help='Celltypes to use')
p$add_argument('--group_variable',    type="character",    help='Group label')
p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per group')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# io$basedir <- file.path(io$basedir,"test")
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
# args$groupA <- "Epiblast"
# args$groupB <- "Paraxial_mesoderm"
# args$group_variable <- "celltype"
# args$celltypes <- "all"
# args$min_cells <- 25
# args$outfile <- file.path(io$basedir,sprintf("results/rna/differential/metacells/celltype/%s_vs_%s.txt.gz",args$groupA,args$groupB))
## END TEST

dir.create(dirname(args$outfile), showWarnings=F, recursive=T)

#####################
## Define settings ##
#####################

# Load utils
source(here::here("rna/differential/utils.R"))

# Define cell types
if (args$celltypes[1]=="all") {
  args$celltypes <- opts$celltypes
} else {
  stopifnot(args$celltypes%in%c(opts$celltypes,"Erythroid","Blood_progenitors"))
}

# Define groups
opts$groups <- c(args$groupA,args$groupB)

opts$min_detection_rate_per_group <- 0.30

# stupid stuff but otherwise the snakemake pipeline doesn't work
if (args$groupA==args$groupB) {
  out <- data.table(gene=NA, logFC=NA, padj_fdr=NA, groupA_N=NA, groupB_N=NA, detection_rate_groupA=NA, detection_rate_groupB=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  stop("groupA and groupB are the same, saving an empty file...")
}

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[celltype%in%args$celltypes]

stopifnot(args$group_variable%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  setnames(args$group_variable,"group") %>%
  # .[,group:=eval(as.name(args$group_variable))] %>%
  .[group%in%opts$groups] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

# Check number of cells per group
tmp <- table(sample_metadata$group)
if (any(tmp<=args$min_cells)) {
  warning("Not enough cells per group to perform DE, saving an empty file...")
  out <- data.table(gene=NA, logFC=NA, padj_fdr=NA, groupA_N=tmp[args$groupA], groupB_N=tmp[args$groupB], detection_rate_groupA=NA, detection_rate_groupB=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  quit(status=0)
} 


#########################
## Load RNA expression ##
#########################

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(
  file = args$sce, 
  normalise = TRUE, 
  cells = sample_metadata$cell
)
sce$group <- sample_metadata$group

stopifnot("logcounts"%in%names(assays(sce)))

#########################################
## Calculate average expression levels ##
#########################################

# calculate detection rate per gene
expr.dt <- data.table(
  gene = rownames(sce),
  mean_groupA = rowMeans(logcounts(sce[,sce$group==args$groupA])) %>% round(2),
  mean_groupB = rowMeans(logcounts(sce[,sce$group==args$groupB])) %>% round(2)
  # detection_rate_groupA = rowMeans(logcounts(sce[,sce$group==args$groupsA])>0) %>% round(2),
  # detection_rate_groupB = rowMeans(logcounts(sce[,sce$group==args$groupsB])>0) %>% round(2)
)

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, opts$groups, opts$min_detection_rate_per_group) %>%
  .[,c("groupA_N","groupB_N"):=list(table(sample_metadata$group)[1],table(sample_metadata$group)[2])]%>% 
  merge(expr.dt, all.y=T, by="gene") %>%
  setorder(padj_fdr, na.last=T)

# Parse columns
out[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3), round(logFC,3))]

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  
