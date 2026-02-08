source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(edgeR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='Cell metadata file')
p$add_argument('--sce',    type="character",    help='SingleCellExperiment file')
p$add_argument('--samples',    type="character",  nargs="+",  help='Samples')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--group_variable',    type="character",    help='Group label')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
io$metadata
args$metadata <- file.path(io$basedir,"processed/rna/sample_metadata.txt")
args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
# args$group_variable <- "Treatment"; args$groupA <- "Control"; args$groupB <- "TG"
args$group_variable <- "Treatment"; args$groupA <- "Control"; args$groupB <- "TG"
# args$samples <- c("p8_Control", "p11_Control", "p8_TG", "p11_TG")
args$samples <- c("p24_Control", "p29_Control", "p24_TG", "p29_TG")
args$outfile <- file.path(io$basedir,sprintf("results/rna/differential/HR3/diff_%s_vs_%s_late_passage.txt.gz",args$groupA,args$groupB))
## END TEST

#####################
## Define settings ##
#####################

# I/O
dir.create(dirname(args$outfile), showWarnings=F, recursive=T)

# Define groups
opts$groups <- c(args$groupA,args$groupB)

if (is.null(args$samples)) {
  args$samples <- opts$groups
}
  
##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(args$metadata) %>%
  .[Group%in%args$samples]

stopifnot(args$group_variable%in%colnames(sample_metadata.dt))

sample_metadata.dt <- sample_metadata.dt %>%
  setnames(args$group_variable,"group") %>%
  # .[,group:=eval(as.name(args$group_variable))] %>%
  .[group%in%opts$groups] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

table(sample_metadata.dt$group)

#########################
## Load RNA expression ##
#########################

# Load SummarizedExperiment object
sce <- readRDS(args$sce)[,sample_metadata.dt$sample]
# rownames(sce) <- rowData(sce)$gene_name
rowData(sce)$gene_name <- rownames(sce)

# colData
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame

# normalise
assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(assay(sce,"counts")),"/"))+1)

#########################################
## Calculate average expression levels ##
#########################################

expr.dt <- data.table(
  gene = rownames(sce),
  transcript_type = rowData(sce)$transcript_type,
  mean_groupA = rowMeans(assay(sce[,sce$group==args$groupA],"logcounts")) %>% round(2),
  mean_groupB = rowMeans(assay(sce[,sce$group==args$groupB],"logcounts")) %>% round(2)
)

#######################
## Feature selection ##
#######################

opts$min.expr <- 5 # 2**5 = 32, at least an average of 32 counts per milion for one group

genes.to.use <- expr.dt[mean_groupA>=opts$min.expr | mean_groupB>=opts$min.expr,gene]

################################################
## Differential expression testing with edgeR ##
################################################

# Convert sce to DGEList
se_edger <- SE2DGEList(sce[genes.to.use,])

# Define design matrix (with intercept)
# library_size <- colSums(assay(sce,"counts")>0)
# design <- model.matrix(~library_size+sce$donor+sce$group)
design <- model.matrix(~sce$group)

# Estimate dispersions
se_edger  <- estimateDisp(se_edger,design)

# Fit GLM
fit <- glmQLFit(se_edger,design)

# Likelihood ratio test
lrt <- glmQLFTest(fit)

# Construct output data.frame
out <- topTags(lrt, n=nrow(lrt))$table %>% as.data.table(keep.rownames=F) %>%
  setnames(c("gene","transcript_type","logFC","logCPM","LR","p.value","padj_fdr")) %>%
  .[,c("logCPM","LR","p.value"):=NULL] %>%
  .[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3), round(logFC,3))] %>%
  .[,c("groupA_N","groupB_N"):=list(table(sample_metadata.dt$group)[1],table(sample_metadata.dt$group)[2])]%>% 
  merge(expr.dt, by=c("gene","transcript_type"), all.y=TRUE) %>%
  setorder(padj_fdr, na.last=T)

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
