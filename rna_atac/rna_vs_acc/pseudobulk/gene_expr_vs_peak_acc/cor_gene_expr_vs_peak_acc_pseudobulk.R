here::i_am("rna_atac/rna_vs_acc/pseudobulk/gene_expr_vs_peak_acc/cor_gene_expr_vs_peak_acc_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
p$add_argument('--peak2gene',  type="character",              help='Peak2gene file') 
p$add_argument('--distance',  type="integer",            default=1e5,      help='Maximum distance for a linkage between a peak and a gene')
p$add_argument('--outfile',          type="character",                help='Output file')
p$add_argument('--test',      action = "store_true",                       help='Testing mode')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped.rds") # io$rna.pseudobulk.sce
# args$atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype.mapped/pseudobulk_PeakMatrix_summarized_experiment.rds")
# args$peak2gene <- file.path(io$basedir,"results/atac/archR/peak_calling/peaks2genes/peaks2genes_all.txt.gz")
# args$distance <- 5e4
# args$outfile <- file.path(io$basedir,"results/rna_atac/gene_expr_vs_peak_acc/pseudobulk/cor_gene_expr_vs_peak_acc_pseudobulk.txt.gz")
# args$test <- TRUE
## END TEST ##

# I/O
dir.create(dirname(args$outfile), showWarnings = F, recursive = TRUE)

##################################
## Load pseudobulk RNA and ATAC ##
##################################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_peakMatrix_pseudobulk.se <- readRDS(args$atac_peak_matrix)

# Normalise ATAC data
assayNames(atac_peakMatrix_pseudobulk.se) <- "counts"
assay(atac_peakMatrix_pseudobulk.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_pseudobulk.se),2,colSums(assay(atac_peakMatrix_pseudobulk.se),na.rm=T),"/"))+1)

# Make sure that samples are consistent
celltypes <- intersect(colnames(rna_pseudobulk.sce),colnames(atac_peakMatrix_pseudobulk.se))
rna_pseudobulk.sce <- rna_pseudobulk.sce[,celltypes]
atac_peakMatrix_pseudobulk.se <- atac_peakMatrix_pseudobulk.se[,celltypes]

#############################
## Load peak2gene linkages ##
#############################

peak2gene.dt <- fread(args$peak2gene) %>%
  .[dist<=args$distance]

# peak2gene.dt <- fread(io$archR.peak2gene.nearest) %>%
#   .[dist<=args$distance]

###########################################################
## Correlate peak accessibility with gene RNA expression ##
###########################################################

# TFs <- intersect(colnames(motifmatcher.se),rownames(rna_pseudobulk.sce))# %>% head(n=5)
# TFs <- c("GATA1","TAL1")
genes <- intersect(rownames(rna_pseudobulk.sce),peak2gene.dt$gene)

if (args$test) {
  genes <- genes %>% head(n=5)
}

cor.dt <- genes %>% map(function(i) {
    tmp <- peak2gene.dt[gene==i,c("gene","peak")]
    print(sprintf("%s (%d/%d)",i,match(i,genes),length(genes)))
    # calculate correlations
    corr_output <- psych::corr.test(
      x = t(logcounts(rna_pseudobulk.sce[i,])), 
      y = t(assay(atac_peakMatrix_pseudobulk.se[tmp$peak,])), 
      ci = FALSE
    )
    
    data.table(
      gene = i,
      peak = tmp$peak,
      cor = corr_output$r[1,] %>% round(3),
      pvalue = corr_output$p[1,] %>% format(digits=3)
    )
}) %>% rbindlist

##########
## Save ##
##########

fwrite(cor.dt, args$outfile, sep="\t")
