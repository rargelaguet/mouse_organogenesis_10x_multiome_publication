source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(edgeR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='Cell metadata file')
p$add_argument('--sce',    type="character",    help='SingleCellExperiment file')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--group_variable',    type="character",    help='Group label')
p$add_argument('--celltypes',    type="character",    help='Celltypes')
# p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per group')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$metadata <- file.path(io$basedir,"results/metacells/per_condition/metacells_metadata.txt.gz")
# args$sce <- file.path(io$basedir,"results/metacells/per_condition/SingleCellExperiment_metacells.rds")
# args$groupA <- "Saline_Control"; args$groupB <- "Bleomycin_day11"; args$group_variable <- "condition"; args$celltypes <- c("fibroblasts")
# args$outfile <- file.path(io$basedir,sprintf("results/differential/metacells/%s/%s_vs_%s.txt.gz",args$group_variable,args$groupA,args$groupB))
## END TEST

opts$groups <- c(args$groupA,args$groupB)

dir.create(dirname(args$outfile), showWarnings=F, recursive=T)

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(args$metadata) %>% setnames("metacell","cell")

stopifnot(args$group_variable%in%colnames(cell_metadata.dt))

cell_metadata.dt <- cell_metadata.dt %>%
  setnames(args$group_variable,"group") %>%
  # .[,group:=eval(as.name(args$group_variable))] %>%
  .[cell_type%in%args$celltypes] %>%
  .[group%in%opts$groups] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

table(cell_metadata.dt$group)

#########################
## Load RNA expression ##
#########################

sce <- load_SingleCellExperiment(args$sce, cells=cell_metadata.dt$cell, normalise = T)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

#########################################
## Calculate average expression levels ##
#########################################

expr.dt <- data.table(
  gene = rownames(sce),
  transcript_type = rowData(sce)$transcript_type,
  mean_groupA = rowMeans(assay(sce[,sce$group==args$groupA],"logcounts")) %>% round(2),
  mean_groupB = rowMeans(assay(sce[,sce$group==args$groupB],"logcounts")) %>% round(2)
) %>% .[,logFC:=mean_groupB-mean_groupA]
  
#######################
## Feature selection ##
#######################

opts$min.expr <- 3 # in logcount space

genes.to.use <- expr.dt[mean_groupA>=opts$min.expr | mean_groupB>=opts$min.expr,gene]

stopifnot(!duplicated(rownames(sce)))

################################################
## Differential expression testing with edgeR ##
################################################

# Convert se to DGEList
se_edger <- SE2DGEList(sce[genes.to.use,])

# Define design matrix (with intercept)
# library_size <- colSums(assay(se,"counts")>0)
# design <- model.matrix(~library_size+se$donor+se$group)
design <- model.matrix(~sce$group)

# Estimate dispersions
se_edger  <- estimateDisp(se_edger,design)

# Fit GLM
fit <- glmQLFit(se_edger,design)

# Likelihood ratio test
lrt <- glmQLFTest(fit)

# Construct output data.frame
out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=T) %>%
  .[,c("rn","FDR")] %>% setnames(c("gene","padj_fdr")) %>%
  merge(expr.dt, by=c("gene"), all.y=TRUE) %>%
  .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3), round(logFC,3))] %>%
  setorder(padj_fdr, na.last=T)

out[is.na(padj_fdr),padj_fdr:=1]

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  
