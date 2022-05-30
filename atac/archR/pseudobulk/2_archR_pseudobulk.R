# https://www.ArchRProject.com/bookdown/how-does-archr-make-pseudo-bulk-replicates.html
here::i_am("atac/archR/pseudobulk/2_archR_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--matrices',     type="character",       nargs="+",   help='Matrices to pseudobulk')
# p$add_argument('--min_cells',     type="integer",    default=50,    help='Minimum number of cells per group')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--outdir',     type="character",    help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$archr_directory <- file.path(io$basedir,"test/processed/atac/archR")
# args$metadata <- file.path(io$basedir,"test/results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$group_by <- "celltype.mapped"
# args$matrices <- "PeakMatrix" # c("PeakMatrix", "GeneScoreMatrix_TSS")
# args$threads <- 1
# # args$min_cells <- 50
# args$outdir <- file.path(io$basedir,sprintf("test/results/atac/archR/pseudobulk/%s/test",args$group_by))
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

########################
## Load cell metadata ##
########################

if (grepl("genotype",args$group_by)) {
  sample_metadata <- fread(args$metadata) %>%
    .[,celltype_genotype:=sprintf("%s-%s",celltype,genotype)] %>%
    .[pass_atacQC==TRUE & doublet_call==FALSE]
} else {
  sample_metadata <- fread(args$metadata) %>%
    .[pass_atacQC==TRUE & doublet_call==FALSE & genotype=="WT"]
}
stopifnot(args$group_by%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$group_by]])]

table(sample_metadata[[args$group_by]])

######################################
## Filter groups by number of cells ##
######################################

# tmp <- table(sample_metadata[[args$group_by]])
# groups.to.use <- tmp[tmp>=args$min_cells] %>% names
# 
# sample_metadata <- sample_metadata[sample_metadata[[args$group_by]]%in%groups.to.use]

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]

# Sanity checks
stopifnot(args$matrices%in%getAvailableMatrices(ArchRProject))

###########################
## Update ArchR metadata ##
###########################

sample_metadata.to.archr <- sample_metadata %>% 
  .[cell%in%rownames(ArchRProject)] %>% setkey(cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(all(rownames(sample_metadata.to.archr) == rownames(getCellColData(ArchRProject))))
ArchRProject <- addCellColData(
  ArchRProject,
  data = sample_metadata.to.archr[[args$group_by]],
  name = args$group_by,
  cells = rownames(sample_metadata.to.archr),
  force = TRUE
)

table(getCellColData(ArchRProject,args$group_by)[[1]])

###################################################
## Pseudobulk into a SummarizedExperiment object ##
###################################################

# if (is.null(args$matrices)) {
# 	args$matrices <- getAvailableMatrices(ArchRProject)
# }

se_list <- list()
for (i in args$matrices) {
  print(sprintf("Calculating pseudobulk matrix for %s",i))

  # summarise
  se_list[[i]] <- getGroupSE(ArchRProject, groupBy = args$group_by, useMatrix = i, divideN = FALSE)
  
  # rename features
  if (grepl("peak",tolower(i),ignore.case=T)) {
    rownames(se_list[[i]]) <- rowData(se_list[[i]]) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
  }
  if (grepl("gene",tolower(i),ignore.case=T)) {
    rownames(se_list[[i]]) <- rowData(se_list[[i]])$name
  }

  # save
  saveRDS(se_list[[i]], file.path(args$outdir,sprintf("pseudobulk_%s_summarized_experiment.rds",i)))
}

# Save stats
to_save.dt <- data.table(table(sample_metadata[[args$group_by]])) %>% 
  setnames(c("group","N")) %>% 
  setorder(group)# %>% .[,included:=group%in%groups.to.use]
fwrite(to_save.dt, file.path(args$outdir,"stats.txt"), sep="\t", quote = F)
