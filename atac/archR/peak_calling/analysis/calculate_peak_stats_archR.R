suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(sparseMatrixStats))


here::i_am("atac/archR/peak_calling/analysis/calculate_peak_stats_archR.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outfile',     type="character",    help='Output file')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--pseudobulk_peak_matrix_file',             type="character",     help='Precomputed peak matrix file (pseudobulk)')
p$add_argument('--ignore_small_celltypes',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))


## START TEST ##
args$metadata <- file.path(io$basedir,"results_new/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
args$celltype_label <- "celltype.predicted"
args$pseudobulk_peak_matrix_file <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn/pseudobulk_PeakMatrix_summarized_experiment.rds")
args$ignore_small_celltypes <- TRUE
args$outfile <- file.path(io$basedir, "results_new/atac/archR/peak_calling/peak_stats.txt.gz")
## END TEST ##

# I/O

# Options

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE]
stopifnot(args$celltype_label%in%colnames(sample_metadata))
sample_metadata <- sample_metadata[!is.na(sample_metadata[[args$celltype_label]])]

# subset celltypes with sufficient number of cells
if (args$ignore_small_celltypes) {
  opts$celltypes <- names(which(table(sample_metadata[[args$celltype_label]])>30))
  sample_metadata <- sample_metadata %>% .[eval(as.name(args$celltype_label))%in%opts$celltypes]
}

########################
## Load ArchR project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

# Subset
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

########################################
## Fetch single-cell ATAC Peak Matrix ##
########################################

print("Fetching single-cell ATAC Peak Matrix...")

atac_peakMatrix.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix", binarize = FALSE)
dim(atac_peakMatrix.se)

# Define peak names
# row.ranges.dt <- rowData(atac_peakMatrix_pseudobulk.se) %>% as.data.table %>% .[,idx:=sprintf("%s_%s_%s",seqnames,start,end)]
row.ranges.dt <- rowRanges(atac_peakMatrix.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,c("chr","start","end")] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
rownames(atac_peakMatrix.se) <- row.ranges.dt$idx

########################################
## Fetch pseudobulk ATAC Peak Matrix ##
########################################

print("Fetching pseudobulk ATAC Peak Matrix...")

atac_peakMatrix_pseudobulk.se <- readRDS(args$pseudobulk_peak_matrix_file)[,opts$celltypes]

# Define peak names
row.ranges.dt <- rowData(atac_peakMatrix_pseudobulk.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,c("chr","start","end")] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]
rownames(atac_peakMatrix_pseudobulk.se) <- row.ranges.dt$idx

###################
## Sanity checks ##
###################

stopifnot(rownames(atac_peakMatrix.se)==rownames(atac_peakMatrix_pseudobulk.se))

##########################
## Calculate peak stats ##
##########################

print("Calculating peak stats...")

peakStats.pseudobulk.dt <- data.table(
  peak = rownames(atac_peakMatrix_pseudobulk.se),
  var = assay(atac_peakMatrix_pseudobulk.se,"PeakMatrix") %>% sparseMatrixStats::rowVars(.) %>% round(5),
  mean = assay(atac_peakMatrix_pseudobulk.se,"PeakMatrix") %>% Matrix::rowMeans(.) %>% round(5)
)

peakStats.singlecell.dt <- data.table(
  peak = rownames(atac_peakMatrix.se),
  var = assay(atac_peakMatrix.se,"PeakMatrix") %>% sparseMatrixStats::rowVars(.) %>% round(5),
  mean = assay(atac_peakMatrix.se,"PeakMatrix") %>% Matrix::rowMeans(.) %>% round(5)
)

peakStats.dt <- merge(peakStats.singlecell.dt, peakStats.pseudobulk.dt, by="peak", suffixes=c("_singlecell","_pseudobulk"))


##########
## Save ##
##########

fwrite(peakStats.dt, args$outfile, sep="\t")
