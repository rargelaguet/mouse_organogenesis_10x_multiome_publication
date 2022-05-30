here::i_am("rna_atac/rna_vs_acc/metacells/gene_expr_vs_promoter_acc/cor_gene_expr_vs_promoter_acc_metacells.R")

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
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/metacells/all_cells/SingleCellExperiment_metacells.rds")
# args$gene_score_matrix <- file.path(io$basedir,"results/atac/archR/metacells/all_cells/GeneScoreMatrix_TSS/GeneScoreMatrix_TSS_summarized_experiment_metacells.rds")
# args$outfile <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/gene_expr_vs_promoter_acc/cor_gene_expr_vs_promoter_acc_metacells.txt.gz")
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
rna_metacells.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_GeneScoreMatrix_metacells.se <- readRDS(args$gene_score_matrix)

# Normalise ATAC data
assayNames(atac_GeneScoreMatrix_metacells.se) <- "counts"
assay(atac_GeneScoreMatrix_metacells.se,"logcounts") <- log(1e6*(sweep(assay(atac_GeneScoreMatrix_metacells.se),2,colSums(assay(atac_GeneScoreMatrix_metacells.se),na.rm=T),"/"))+1)

# hist(assay(atac_GeneScoreMatrix_metacells.se,"logcounts")[1:1000,])

###########################################
## Convert to long data.tables and merge ##
###########################################

rna_metacells.dt <- logcounts(rna_metacells.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

atac_gene_scores_metacells.dt <- as.matrix(assay(atac_GeneScoreMatrix_metacells.se,"logcounts")) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="gene", value.name="acc")


# Merge
rna_atac.dt <- merge(rna_metacells.dt, atac_gene_scores_metacells.dt, by = c("gene","celltype"))

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
