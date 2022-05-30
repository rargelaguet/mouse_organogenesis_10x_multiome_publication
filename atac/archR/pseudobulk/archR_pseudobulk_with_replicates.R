here::i_am("atac/archR/pseudobulk/archR_pseudobulk_with_replicates.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atac_matrix_file',    type="character",    help='ATAC matrix')
p$add_argument('--atac_matrix_name',    type="character",    help='ATAC matrix')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--nrep',       type="integer",       default=5,      help='Number of replicates per group (cells sampled with replacement)')
p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per replicate')
p$add_argument('--fraction_cells_per_replicate',       type="double",       default=0.3,      help='Percentage of cells per replicate')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--outdir',     type="character",    help='Output directory')
# p$add_argument('--test_mode',       action="store_true",  help='Test mode? subset data')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args$group_by <- "celltype"
# args$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$atac_matrix_name <- "GeneScoreMatrix_TSS"
# args$atac_matrix_file <- file.path(io$basedir,sprintf("processed/atac/archR/Matrices/%s_summarized_experiment.rds",args$atac_matrix_name))
# args$nrep <- 5
# args$min_cells <- 25
# args$fraction_cells_per_replicate <- 0.30
# args$outdir <- file.path(io$basedir,sprintf("results/atac/archR/pseudobulk/%s/%s",args$group_by,args$atac_matrix_name))
# # args$test_mode <- TRUE
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(args$metadata) %>%
  .[,celltype_genotype:=sprintf("%s-%s",celltype,genotype)] %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$group_by)))] %>%
  setnames(args$group_by,"group")

print(cell_metadata.dt[,.N,by="group"])

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

######################
## Load ATAC matrix ##
######################

print("Loading ATAC matrix...")

atac.se <- readRDS(args$atac_matrix_file)[,cell2replicate.dt$cell] %>% as(.,"SingleCellExperiment")
atac.se$group_by <- cell2replicate.dt$replicate
assayNames(atac.se) <- "counts"

# if (args$test_mode) {
#   print("testing mode activated, subsetting number of features...")
#   atac.se <- atac.se[1:100,]
# }

################
## Pseudobulk ##
################

print("Pseudobulking...")

atac_pseudobulk.se <- pseudobulk_sce_fn(
  x = atac.se,
  assay = "counts",
  by = "group_by",
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(atac_pseudobulk.se) <- "counts"

###################
## Normalisation ##
###################

logcounts(atac_pseudobulk.se) <- log(1e6*(sweep(counts(atac_pseudobulk.se),2,colSums(counts(atac_pseudobulk.se)),"/"))+1)

##########
## Save ##
##########

# edit rowData
# if (args$atac_matrix_name=="PeakMatrix") {
#   tmp <- rownames(atac_pseudobulk.se) %>% strsplit(":") %>% map_chr(2)
#   rowData(atac_pseudobulk.se)$start <- tmp %>% strsplit("-") %>% map_chr(1)
#   rowData(atac_pseudobulk.se)$end <- tmp %>% strsplit("-") %>% map_chr(2)
# }


# Edit colData
if (args$group_by=="celltype_genotype") {
  atac_pseudobulk.se$celltype <- colnames(atac_pseudobulk.se) %>% strsplit("-") %>% map_chr(1)
  atac_pseudobulk.se$genotype <- colnames(atac_pseudobulk.se) %>% strsplit("-") %>% map_chr(2) %>% strsplit("_rep") %>% map_chr(1)
  atac_pseudobulk.se$celltype_genotype <- sprintf("%s_%s",atac_pseudobulk.se$celltype,atac_pseudobulk.se$genotype)
} else if (args$group_by=="celltype") {
  atac_pseudobulk.se$celltype <- colnames(atac_pseudobulk.se) %>% strsplit("_rep") %>% map_chr(1)
}

# Save Summarizedexperiment
saveRDS(atac_pseudobulk.se, file.path(args$outdir,sprintf("%s_pseudobulk_with_replicates.rds",args$atac_matrix_name)))

# Save cell2replicate assignment
fwrite(cell2replicate.dt, file.path(args$outdir,"cell2replicate.txt.gz"), sep="\t", quote = F)

# Save stats
# to_save.dt <- cell2replicate.dt[,.N,by="replicate"] %>% setnames(c("group","N")) %>% setorder(group)
# fwrite(to_save.dt, file.path(args$outdir,"stats.txt"), sep="\t", quote = F)

# Create a completion token
# file.create(file.path(args$outdir,"completed.txt"))
