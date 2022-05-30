here::i_am("rna/pseudobulk/create_pseudobulk_metadata_with_replicates.R")


source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--group_by',    type="character",    help='')
p$add_argument('--nrep',       type="integer",       default=5,      help='Number of replicates per group (cells sampled with replacement)')
p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per replicate')
p$add_argument('--percentage_cells_per_replicate',       type="double",       default=0.3,      help='Percentage of cells per replicate')
p$add_argument('--outdir',      type="character",    help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
io$basedir <- file.path(io$basedir,"test")
args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
args$group_by <- "celltype"
args$nrep <- 5
args$min_cells <- 25
args$percentage_cells_per_replicate <- 0.30
args$outfile <- file.path(io$basedir,sprintf("results/rna/pseudobulk/%s/cell2replicate.txt.gz",args$group_by))
## END TEST ##

dir.create(dirname(args$outfile), showWarnings = F, recursive = T)

###################
## Load metadata ##
###################

# Load cell metadata
cell_metadata.dt <- fread(args$metadata) %>%
  .[,celltype_genotype:=sprintf("%s-%s",celltype,genotype)] %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(eval(as.name(args$group_by)))] %>%
  setnames(args$group_by,"group")

print(table(cell_metadata.dt$group))

##################################
## Create pseudobulk replicates ##
##################################

cell2group.dt <- unique(cell_metadata.dt$group) %>% map(function(i) {
  tmp <- cell_metadata.dt[group==i]
  if ((args$percentage_cells_per_replicate*nrow(tmp))<=args$min_cells) {
    ncells_per_replicate <- args$min_cells
  } else {
    ncells_per_replicate <- round(args$percentage_cells_per_replicate*nrow(tmp))
  }
  seq(1,args$nrep) %>% map(function(j) {
    tmp[sample.int(nrow(tmp),ncells_per_replicate)] %>% 
      .[,replicate:=sprintf("%s_rep%s",i,j)] %>%
      .[,c("cell","group","replicate")] %>% 
      return
  }) %>% rbindlist %>% return
}) %>% rbindlist


stats.dt <- cell2group.dt[,.(ncells=.N),c("group","replicate")]
print(stats.dt)

##########
## Save ##
##########

fwrite(cell2group.dt, args$outfile, sep="\t", quote = F)

