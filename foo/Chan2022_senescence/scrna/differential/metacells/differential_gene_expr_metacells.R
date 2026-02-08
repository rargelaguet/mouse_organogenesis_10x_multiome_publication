# here::i_am("rna/differential/differential_gene_expr.R")

# Load default settings
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
# p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per group')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$metadata <- file.path(io$basedir,"results/rna/metacells/scTE/metacells_metadata.txt.gz")
# args$sce <- file.path(io$basedir,"results/rna/metacells/scTE/SingleCellExperiment_metacells.rds")
# args$groupA <- "D16_TRA-1-60_positive"
# args$groupB <- "D16_TRA-1-60_negative"
# args$group_variable <- "alias"
# args$outfile <- file.path(io$basedir,sprintf("results/rna/differential/metacells/diff_%s_%s_vs_%s.txt.gz",args$group_variable,args$groupA,args$groupB))
## END TEST

opts$groups <- c(args$groupA,args$groupB)

dir.create(dirname(args$outfile), showWarnings=F, recursive=T)

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(args$metadata) %>% setnames("metacell","cell")

stopifnot(args$group_variable%in%colnames(cell_metadata.dt))

# stopifnot(!"group"%in%colnames(cell_metadata.dt))

if ("group"%in%colnames(cell_metadata.dt) & args$group_variable!="group") {
  cell_metadata.dt <- cell_metadata.dt %>%  .[,group:=NULL]
}

cell_metadata.dt <- cell_metadata.dt %>% setnames(args$group_variable,"group")

cell_metadata.dt <- cell_metadata.dt %>% 
    .[group%in%opts$groups] %>%
    .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

table(cell_metadata.dt$group)

#########################
## Load RNA expression ##
#########################

# Load SummarizedExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=cell_metadata.dt$cell, normalise = T)
# sce <- sce[!duplicated(rownames(sce))]

# colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

# normalise (cpm)
# assay(se,"logcounts") <- log2(1e6*(sweep(assay(se,"counts"),2,colSums(assay(se,"counts")),"/"))+1)
# assay(sce,"logcounts") <- 1e6*(sweep(assay(sce,"counts"),2,colSums(assay(sce,"counts")),"/"))

#########################################
## Calculate average expression levels ##
#########################################

# calculate detection rate per gene
expr.dt <- data.table(
  gene = rownames(sce),
  transcript_type = rowData(sce)$transcript_type,
  mean_groupA = rowMeans(assay(sce[,sce$group==args$groupA],"logcounts")) %>% round(2),
  mean_groupB = rowMeans(assay(sce[,sce$group==args$groupB],"logcounts")) %>% round(2)
) %>% .[,logFC:=mean_groupB-mean_groupA]
  
#######################
## Feature selection ##
#######################

opts$min.expr <- 2 # in logcount space

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
  setnames(c("gene","transcript_type","logFC","logCPM","LR","p.value","padj_fdr")) %>%
  .[,c("logCPM","LR","p.value","logFC"):=NULL] %>%
  # .[,c("groupA_N","groupB_N"):=list(table(cell_metadata.dt$group)[1],table(cell_metadata.dt$group)[2])]%>% 
  merge(expr.dt, by=c("gene","transcript_type"), all.y=TRUE) %>%
  .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3), round(logFC,3))] %>%
  setorder(padj_fdr, na.last=T)

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  
