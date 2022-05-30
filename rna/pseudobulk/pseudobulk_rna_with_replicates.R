here::i_am("rna/pseudobulk/pseudobulk_rna_with_replicates.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

# suppressPackageStartupMessages(library(Seurat))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',         type="character",    help='SingleCellExperiment object')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--cell2replicate',    type="character",    help='')
p$add_argument('--group_by',    type="character",    help='Metadata column to group cells by')
p$add_argument('--nrep',       type="integer",       default=5,      help='Number of replicates per group (cells sampled with replacement)')
p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per replicate')
p$add_argument('--fraction_cells_per_replicate',       type="double",       default=0.3,      help='Percentage of cells per replicate')
p$add_argument('--outdir',      type="character",    help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
# args$group_by <- "celltype_genotype"
# args$nrep <- 5
# args$min_cells <- 25
# args$fraction_cells_per_replicate <- 0.30
# args$outdir <- file.path(io$basedir,sprintf("results/rna/pseudobulk/%s",args$group_by))
## END TEST ##

dir.create(args$outdir, showWarnings = F, recursive = T)

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(args$metadata) %>%
  .[,celltype_genotype:=as.character(NA)] %>% .[!is.na(genotype),celltype_genotype:=sprintf("%s-%s",celltype,genotype)] %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$group_by)))] %>%
  setnames(args$group_by,"group")

# if (grep("genotype",args$group_by)) {
#   cell_metadata.dt <- cell_metadata.dt[sample%in%c("E8.5_CRISPR_T_KO", "E8.5_CRISPR_T_WT")]
# }

print(table(cell_metadata.dt$group))

# Remove groups with not enough cells
groups.to.remove <- cell_metadata.dt[,.N,by="group"] %>% .[N<=args$min_cells,group]
print(sprintf("Removing the following groups because they have less than %d cells:",args$min_cells))
print(cell_metadata.dt[group%in%groups.to.remove,.N,by="group"])
cell_metadata.dt <- cell_metadata.dt[!group%in%groups.to.remove]

##################################
## Create pseudobulk replicates ##
##################################

cell2replicate.dt <- unique(cell_metadata.dt$group) %>% map(function(i) {
  tmp <- cell_metadata.dt[group==i]
  if ((args$fraction_cells_per_replicate*nrow(tmp))<=args$min_cells) {
    ncells_per_replicate <- args$min_cells
  } else {
    ncells_per_replicate <- round(args$fraction_cells_per_replicate*nrow(tmp))
  }
  seq(1,args$nrep) %>% map(function(j) {
    tmp[sample.int(nrow(tmp),ncells_per_replicate)] %>% 
      .[,replicate:=sprintf("%s_rep%s",i,j)] %>%
      .[,c("cell","group","replicate")] %>% 
      return
  }) %>% rbindlist %>% return
}) %>% rbindlist

stats.dt <- cell2replicate.dt[,.(ncells=.N),c("group","replicate")]
print(stats.dt)

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(args$sce, cells=cell2replicate.dt$cell)
sce$group_by <- cell2replicate.dt$replicate

################
## Pseudobulk ##
################

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = "group_by",
  fun = "sum",
  scale = FALSE
)

assayNames(sce_pseudobulk) <- "counts"

###################
## Normalisation ##
###################

logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)

##########
## Save ##
##########

# Save cell2replicate assignment
fwrite(cell2replicate.dt, file.path(args$outdir,"cell2replicate.txt.gz"), sep="\t", quote = F)

# Save SingleCellExperiment
if (args$group_by=="celltype") {
  sce_pseudobulk$celltype <- colnames(sce_pseudobulk) %>% strsplit("_rep") %>% map_chr(1)
}
saveRDS(sce_pseudobulk, file.path(args$outdir,"SingleCellExperiment_pseudobulk_with_replicates.rds"))

##########
## TEST ##
##########

# sort(logcounts(sce_pseudobulk)["Tfap2c",])
