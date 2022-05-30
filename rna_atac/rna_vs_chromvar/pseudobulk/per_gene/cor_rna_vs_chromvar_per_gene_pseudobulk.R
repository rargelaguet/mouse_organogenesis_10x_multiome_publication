here::i_am("rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromvar_per_gene_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',                type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--atac_chromvar_matrix',  type="character",              help='ATAC chromVAR matrix (pseudobulk)') 
p$add_argument('--motif_annotation',   type="character",              help='Motif annotation') 
p$add_argument('--motif2gene',   type="character",              help='') 
p$add_argument('--outfile',            type="character",              help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$motif_annotation <- "JASPAR"
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
# args$atac_chromvar_matrix <- file.path(io$basedir,sprintf("results/atac/archR/chromvar/pseudobulk/celltype/chromVAR_deviations_%s_pseudobulk_archr.rds",args$motif_annotation))
# args$motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s_motif2gene.txt.gz",args$motif_annotation))
# args$outfile <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/%s/cor_rna_vs_chromvar_%s_per_gene_pseudobulk.txt.gz",args$motif_annotation,args$motif_annotation))
## END TEST ##

# I/O
dir.create(dirname(args$outfile), showWarnings=F, recursive = T)

# Options
# opts$motif_annotation <- args$motif_annotation

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_chromvar_pseudobulk.se <- readRDS(args$atac_chromvar_matrix)

################
## Parse data ##
################

# Load motif2gene annotation
motif2gene.dt <- fread(args$motif2gene) %>%
  .[motif%in%rownames(atac_chromvar_pseudobulk.se) & gene%in%toupper(rownames(rna_pseudobulk.sce))]# %>%
  # .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

# Duplicated motifs
print(sprintf("%s genes with duplicated motifs.",nrow(motif2gene.dt[,.N,c("gene","motif")][N>1])))

# Subset TFs
rna_tf_pseudobulk.sce <- rna_pseudobulk.sce[str_to_title(motif2gene.dt$gene),]
rownames(rna_tf_pseudobulk.sce) <- toupper(rownames(rna_tf_pseudobulk.sce))
atac_chromvar_pseudobulk.se <- atac_chromvar_pseudobulk.se[motif2gene.dt$motif,]

###########################################
## Convert to long data.tables and merge ##
###########################################

rna_tf_pseudobulk.dt <- logcounts(rna_tf_pseudobulk.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

atac_chromvar_pseudobulk.dt <- as.matrix(assay(atac_chromvar_pseudobulk.se,"z")) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="motif", value.name="chromvar") %>%
  merge(motif2gene.dt[,c("motif","gene")], by="motif")

# Merge
rna_atac.dt <- merge(rna_tf_pseudobulk.dt, atac_chromvar_pseudobulk.dt, by = c("gene","celltype"))

# Print stats
print(sprintf("Number of TFs: %s",length(unique(rna_atac.dt$gene))))
print(sprintf("Number of celltypes: %s",length(unique(rna_atac.dt$celltype))))

##########################
## Correlation analysis ##
##########################

# opts$threshold_fdr <- 0.10

cor.dt <- rna_atac.dt %>% copy %>%
  .[,c("chromvar","expr"):=list(chromvar + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[, .(V1 = unlist(cor.test(chromvar, expr)[c("estimate", "p.value")])), by = c("gene","motif")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene+motif ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  # .[, sig := p<=opts$threshold_fdr] %>% 
  setorder(p, na.last = T)

# Save
fwrite(cor.dt, args$outfile, sep="\t", quote=F)
