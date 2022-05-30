here::i_am("atac/archR/feature_stats/archR_calculate_feature_stats.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(sparseMatrixStats))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outfile',     type="character",    help='Output file')
p$add_argument('--cells_metadata', type="character", help='Cell metadata')
p$add_argument('--metacells_metadata', type="character", help='Metacell metadata')
# p$add_argument('--pseudobulk_metadata', type="character", help='Pseudobulk metadata')
p$add_argument('--matrix',             type="character",     help='Matrix name')
p$add_argument('--group_by',             type="character",     help='Metadata column to group by')
p$add_argument('--atac_matrix_cells', type="character", help='ATAC matrix (cells)')
p$add_argument('--atac_matrix_metacells', type="character", help='ATAC matrix (metacells)')
p$add_argument('--atac_matrix_pseudobulk', type="character", help='ATAC matrix (pseudobulk)')
p$add_argument('--ignore_small_groups',  action="store_true",  help='')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args$cells_metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$metacells_metadata <- file.path(io$basedir,"results/atac/archR/metacells/GeneScoreMatrix_TSS/metacells_metadata.txt.gz")
# # args$pseudobulk_metadata <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype/stats.txt")
# args$matrix <- "GeneScoreMatrix_TSS"
# args$group_by <- "celltype"
# args$atac_matrix_cells <- file.path(io$basedir,sprintf("processed/atac/archR/Matrices/%s_summarized_experiment.rds",args$matrix))
# args$atac_matrix_metacells <- file.path(io$basedir,sprintf("results/atac/archR/metacells/%s/%s_summarized_experiment_metacells.rds",args$matrix,args$matrix))
# args$atac_matrix_pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/pseudobulk/celltype/%s/pseudobulk_%s_summarized_experiment.rds",args$matrix,args$matrix))
# args$ignore_small_groups <- TRUE
# args$outfile <- file.path(io$basedir, sprintf("results/atac/archR/feature_stats/%s/%s_%s_feature_stats.txt.gz",args$matrix,args$matrix,args$group_by))
## END TEST ##

dir.create(dirname(args$outfile), showWarnings=F)

########################
## Load cell metadata ##
########################

if (grepl("genotype",args$group_by)) {
  cell_metadata.dt <- fread(args$cells_metadata) %>%
    .[,celltype_genotype:=sprintf("%s-%s",celltype,genotype)] %>%
    .[pass_atacQC==TRUE & doublet_call==FALSE]
} else {
  cell_metadata.dt <- fread(args$cells_metadata) %>%
    .[pass_atacQC==TRUE & doublet_call==FALSE & genotype=="WT"]
}

stopifnot(args$group_by%in%colnames(cell_metadata.dt))
cell_metadata.dt <- cell_metadata.dt[!is.na(cell_metadata.dt[[args$group_by]])]

############################
## Load metacell metadata ##
############################

if (grepl("genotype",args$group_by)) {
  metacell_metadata.dt <- fread(args$metacells_metadata) %>%
    .[,celltype_genotype:=sprintf("%s-%s",celltype,genotype)]
} else {
  metacell_metadata.dt <- fread(args$metacells_metadata) %>% .[genotype=="WT"]
}

stopifnot(args$group_by%in%colnames(metacell_metadata.dt))
metacell_metadata.dt <- metacell_metadata.dt[!is.na(metacell_metadata.dt[[args$group_by]])]

##############################
## Load pseudobulk metadata ##
##############################

# pseudobulk_metadata.dt <- fread(args$pseudobulk_metadata)

###################
## Filter groups ##
###################

group_stats.dt <- merge(
  cell_metadata.dt[,.(n_cells=.N),by=eval(args$group_by)],
  metacell_metadata.dt[,.(n_metacells=.N),by=eval(args$group_by)],
  by = args$group_by
)

# subset groups with sufficient number of cells
opts$min_cells <- 50
opts$min_metacells <- 5
if (args$ignore_small_groups) {
  groups.to.use <- group_stats.dt[n_cells>=opts$min_cells & n_metacells>=opts$min_metacells][[args$group_by]]
  print(sprintf("Groups removed because they contain less than %s cells and %s metacells: %s",opts$min_cells,opts$min_metacells, paste(setdiff(unique(cell_metadata.dt[[args$group_by]]),groups.to.use), collapse=" ") ))
  cell_metadata.dt <- cell_metadata.dt %>% .[eval(as.name(args$group_by))%in%groups.to.use]
  metacell_metadata.dt <- metacell_metadata.dt %>% .[eval(as.name(args$group_by))%in%groups.to.use]
}

######################
## Load ATAC Matrix ##
######################

# cells
print(sprintf("Fetching single-cell ATAC %s matrix...",args$matrix))
atac_cells.se <- readRDS(args$atac_matrix_cells)[,cell_metadata.dt$cell]
atac_cells.se

# metacells
print(sprintf("Fetching metacells ATAC %s matrix...",args$matrix))
atac_metacells.se <- readRDS(args$atac_matrix_metacells)[,metacell_metadata.dt$metacell]
atac_metacells.se

# pseudobulk
print(sprintf("Fetching pseudobulk ATAC %s matrix...",args$matrix))
atac_pseudobulk.se <- readRDS(args$atac_matrix_pseudobulk)[,groups.to.use]
atac_pseudobulk.se

# Sanity checks
stopifnot(rownames(atac_cells.se)==rownames(atac_pseudobulk.se))
stopifnot(rownames(atac_cells.se)==rownames(atac_metacells.se))

####################
## Normalise ATAC ##
####################

print("Normalising ATAC data...")
assay(atac_metacells.se) <- 1e6*(sweep(assay(atac_metacells.se),2,colSums(assay(atac_metacells.se),na.rm=T),"/"))
assay(atac_metacells.se) <- log2(assay(atac_metacells.se)+0.5)

assay(atac_pseudobulk.se) <- 1e6*(sweep(assay(atac_pseudobulk.se),2,colSums(assay(atac_pseudobulk.se),na.rm=T),"/"))
assay(atac_pseudobulk.se) <- log2(assay(atac_pseudobulk.se)+0.5)

##########################
## Calculate peak stats ##
##########################

print("Calculating feature stats...")

feature_stats_cells.dt <- data.table(
  feature = rownames(atac_cells.se),
  var_cells = sparseMatrixStats::rowVars(assay(atac_cells.se)) %>% round(5),
  mean_cells = Matrix::rowMeans(assay(atac_cells.se)) %>% round(5)
  # min_cells = apply(assay(atac_cells.se),1,min),
  # max_cells = apply(assay(atac_cells.se),1,max)
) %>% setkey(feature)

feature_stats_metacells.dt <- data.table(
  feature = rownames(atac_metacells.se),
  var_metacells = apply(assay(atac_metacells.se),1,var) %>% round(5),
  mean_metacells = apply(assay(atac_metacells.se),1,mean) %>% round(5)
  # min_metacells = apply(assay(atac_metacells.se),1,min) %>% round(5),
  # max_metacells = apply(assay(atac_metacells.se),1,max) %>% round(5)
) %>% setkey(feature)

feature_stats_pseudobulk.dt <- data.table(
  feature = rownames(atac_pseudobulk.se),
  var_pseudobulk = apply(assay(atac_pseudobulk.se),1,var) %>% round(5),
  mean_pseudobulk = apply(assay(atac_pseudobulk.se),1,mean) %>% round(5)
  # min_pseudobulk = apply(assay(atac_pseudobulk.se),1,min) %>% round(5),
  # max_pseudobulk = apply(assay(atac_pseudobulk.se),1,max) %>% round(5)
) %>% setkey(feature)

# feature_stats.dt <- merge(feature_stats_cells.dt, feature_stats_pseudobulk.dt, by="feature", suffixes=c("_singlecell","_pseudobulk"))
feature_stats.dt <- Reduce(data.table::merge.data.table,list(feature_stats_cells.dt, feature_stats_metacells.dt, feature_stats_pseudobulk.dt))

print(head(feature_stats.dt))

##########
## Save ##
##########

fwrite(feature_stats.dt, args$outfile, sep="\t")

##########
## TEST ##
##########

# to.plot <- feature_stats.dt %>% head(n=1e4)
# to.plot <- feature_stats.dt[var_pseudobulk>=0.10 & var_metacells<1]
# 
# ggscatter(to.plot, x="var_metacells", y="var_pseudobulk", size=1) + 
#   geom_abline(slope=1, intercept=0)
# 
# ggscatter(to.plot, x="var_cells", y="var_pseudobulk", size=1) + 
#   geom_abline(slope=1, intercept=0)
# 
# 
# 
# ggscatter(to.plot, x="mean_metacells", y="var_metacells", size=1) + 
#   stat_smooth(method="loess")
# 
# ggscatter(to.plot, x="mean_pseudobulk", y="var_pseudobulk", size=0.5) + 
#   stat_smooth(method="loess")
