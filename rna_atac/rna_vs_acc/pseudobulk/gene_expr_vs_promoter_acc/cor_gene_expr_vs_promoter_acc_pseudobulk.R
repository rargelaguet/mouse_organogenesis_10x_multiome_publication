here::i_am("rna_atac/rna_vs_acc/pseudobulk/gene_expr_vs_promoter_acc/cor_gene_expr_vs_promoter_acc_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--gene_score_matrix',  type="character",              help='ATAC Gene score matrix (pseudobulk)') 
p$add_argument('--outfile',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped.rds") 
# args$gene_score_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype.mapped/pseudobulk_GeneScoreMatrix_TSS_summarized_experiment.rds") # io$archR.pseudobulk.GeneMatrix.se
# args$outfile <- file.path(io$basedir,"results/rna_atac/gene_expr_vs_promoter_acc/pseudobulk/cor_gene_expr_vs_promoter_acc_pseudobulk.txt.gz")
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
dir.create(dirname(args$outfile), showWarnings=FALSE, recursive=TRUE)

# Options

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_pseudobulk_GeneScoreMatrix.se <- readRDS(args$gene_score_matrix)
assayNames(atac_pseudobulk_GeneScoreMatrix.se) <- "counts"

# Normalise ATAC data
assay(atac_pseudobulk_GeneScoreMatrix.se,"logcounts") <- log(1e6*(sweep(assay(atac_pseudobulk_GeneScoreMatrix.se),2,colSums(assay(atac_pseudobulk_GeneScoreMatrix.se),na.rm=T),"/"))+1)

###########################################
## Convert to long data.tables and merge ##
###########################################

rna_pseudobulk.dt <- logcounts(rna_pseudobulk.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

atac_gene_scores_pseudobulk.dt <- as.matrix(assay(atac_pseudobulk_GeneScoreMatrix.se)) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="gene", value.name="acc")

# Merge
rna_atac.dt <- merge(rna_pseudobulk.dt, atac_gene_scores_pseudobulk.dt, by = c("gene","celltype"))

##########################
## Correlation analysis ##
##########################

cor.dt <- rna_atac.dt %>% copy %>%
  .[,c("acc","expr"):=list(acc + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>%
  .[, .(V1 = unlist(cor.test(acc, expr)[c("estimate", "p.value")])), by = c("gene")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  # .[, sig := padj_fdr<=0.10] %>% 
  setorder(padj_fdr, na.last = T)

cor.dt[,c("p","r","padj_fdr"):=list(format(p,digits=3),round(r,3), format(padj_fdr,digits=3))]

# Save
fwrite(cor.dt, args$outfile, sep="\t", quote=F)
