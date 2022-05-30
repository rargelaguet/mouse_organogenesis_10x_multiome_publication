here::i_am("rna/differential/metacells/differential.R")

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
p$add_argument('--test_mode', action="store_true", help='Test mode? subset number of cells')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$metadata <- file.path(io$basedir,"results/rna/metacells/all_cells/metacells_metadata.txt.gz")
# args$sce <- file.path(io$basedir,"results/rna/metacells/all_cells/SingleCellExperiment_metacells.rds")
# # args$groupA <- "Epiblast"; args$groupB <- "Paraxial_mesoderm"; args$group_variable <- "celltype"; args$celltypes <- "all"
# args$groupA <- "WT"; args$groupB <- "T_KO"; args$group_variable <- "genotype"; args$celltypes <- "ExE_endoderm"
# args$test_mode <- FALSE
# args$min_cells <- 5
# args$outfile <- NULL
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

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(args$metadata) %>%
  .[celltype%in%args$celltypes]

stopifnot(args$group_variable%in%colnames(cell_metadata.dt))

cell_metadata.dt <- cell_metadata.dt %>%
  setnames(args$group_variable,"group") %>%
  # .[,group:=eval(as.name(args$group_variable))] %>%
  .[group%in%opts$groups] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

if (isTRUE(args$test_mode)) {
  print("Testing mode activated")
  cell_metadata.dt <- cell_metadata.dt %>% split(.,.$group) %>% map(~ head(.,n=250)) %>% rbindlist
}

table(cell_metadata.dt$group)

# Filter groups with small number of cells
tmp <- table(cell_metadata.dt$group)
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
  cells = cell_metadata.dt$metacell
)
sce$group <- cell_metadata.dt$group

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

#######################
## Feature selection ##
#######################

opts$min.expr <- 3 # 2**3 = 8, at least an average of 4 counts per milion for each group

genes.to.use <- expr.dt[mean_groupA>=opts$min.expr | mean_groupB>=opts$min.expr,gene]

################################################
## Differential expression testing with edgeR ##
################################################

# Convert SCE to DGEList
sce_edger <- scran::convertTo(sce[genes.to.use,], type="edgeR")

# Define design matrix (with intercept)
cdr <- colSums(counts(sce)>0)
design <- model.matrix(~cdr+sce$group)

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
  .[,c("groupA_N","groupB_N"):=list(table(cell_metadata.dt$group)[1],table(cell_metadata.dt$group)[2])]%>% 
  merge(expr.dt, by="gene", all.y=TRUE) %>%
  setorder(padj_fdr, na.last=T)

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  
