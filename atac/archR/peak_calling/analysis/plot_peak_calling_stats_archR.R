# TO-DO: INCORPORATE METACELLS

here::i_am("atac/archR/peak_calling/analysis/plot_peak_calling_stats_archR.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(sparseMatrixStats))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--peak_matrix_file',             type="character",     help='Precomputed peak matrix file')
p$add_argument('--min_peak_score',     type="integer",    default=10,    help='Minimum peak score')
p$add_argument('--pseudobulk_peak_matrix_file',             type="character",     help='Precomputed peak matrix file (pseudobulk)')
p$add_argument('--atac_peak_metadata',             type="character",     help='Peak metadata')
p$add_argument('--atac_peak_stats',             type="character",     help='Peak stats')
p$add_argument('--outdir',     type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# args$peak_matrix_file <- file.path(io$basedir,"processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
# args$pseudobulk_peak_matrix_file <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype/PeakMatrix/pseudobulk_PeakMatrix_summarized_experiment.rds")
# args$min_peak_score <- 10
# args$atac_peak_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
# args$atac_peak_stats <- file.path(io$basedir,"results/atac/archR/feature_stats/PeakMatrix/PeakMatrix_celltype_stats.txt.gz")
# args$outdir <- file.path(io$basedir, "results/atac/archR/peak_calling")
## END TEST ##

#####################
## Define settings ##
#####################

# I/O

# Options

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  # .[pass_atacQC==TRUE & doublet_call==FALSE]
  .[pass_atacQC==TRUE]

##########################
## Load peak annotation ##
##########################

peak_metadata.dt <- fread(args$atac_peak_metadata) %>%
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,peakType:=factor(peakType,levels=c("Promoter","Intronic","Exonic","Distal"))]

# Load peak stats
peak_stats.dt <- fread(args$atac_peak_stats) %>% setnames("feature","peak")

stopifnot(c("mean_cells","var_cells","mean_pseudobulk","var_pseudobulk")%in%colnames(peak_stats.dt))

########################################
## Fetch single-cell ATAC Peak Matrix ##
########################################

print("Loading single-cell ATAC Peak Matrix...")

atac_peakMatrix_cells.se <- readRDS(args$peak_matrix_file)[,sample_metadata$cell]

########################################
## Fetch pseudobulk ATAC Peak Matrix ##
########################################

print("Fetching pseudobulk ATAC Peak Matrix...")

atac_peakMatrix_pseudobulk.se <- readRDS(args$pseudobulk_peak_matrix_file)[,opts$celltypes]

# Normalise
# assay(atac_peakMatrix_pseudobulk.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_pseudobulk.se),2,colSums(assay(atac_peakMatrix_pseudobulk.se),na.rm=T),"/"))+1)

###################
## Sanity checks ##
###################

stopifnot(rownames(atac_peakMatrix_cells.se)==rownames(atac_peakMatrix_pseudobulk.se))

####################################################
## Plot distribution of ATAC fragments per sample ##
####################################################

to.plot <- unique(sample_metadata$sample) %>% map(function(i) {
  mtx <- as.matrix(assay(atac_peakMatrix_cells.se[,sample_metadata[sample==i,cell]])[,1:500])
  data.table(table(mtx[mtx<=10])) %>% setnames(c("value","N")) %>% .[,sample:=i]
}) %>% rbindlist

to.plot[N>=1e7,N:=1e7]

p <- ggbarplot(to.plot, x="value", y="N", fill="gray70") +
  facet_wrap(~sample) +
  theme(
    axis.text =  element_text(size=rel(0.65)),
    axis.title.x = element_blank()
  )

pdf(file.path(args$outdir,"atac_fragment_distribution_per_sample.pdf"), width=8, height=5)
print(p)
dev.off()

###############################
## Plot mean versus variance ##
###############################

to.plot <- peak_metadata.dt[chr=="chr1",c("peak","peakType","score")] %>%
  .[score>=args$min_peak_score] %>%
  merge(peak_stats.dt,by="peak")

p1 <- ggscatter(to.plot, x="mean_cells", y="var_cells", size=1) +
  labs(x="Mean (single-cell)", y="Variance (single-cell)") +
  theme(
    axis.text = element_text(size=rel(0.7))
  )

p2 <- ggscatter(to.plot, x="mean_pseudobulk", y="var_pseudobulk", size=1) +
  labs(x="Mean (pseudobulk)", y="Variance (pseudobulk)") +
  theme(
    axis.text = element_text(size=rel(0.7))
  )

p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)

pdf(file.path(args$outdir,"atac_peak_mean_vs_variance.pdf"), width = 8, height = 4.5)
print(p)
dev.off()

####################################################
## Plot minimum peak score versus number of peaks ##
####################################################

to.plot <- seq(1,250,by=1) %>% map( function(x) {
    data.table(N=peak_metadata.dt[score>=x,.N], min_score=x)
  }) %>% rbindlist

p <- ggscatter(to.plot, x="min_score", y="N", size=0.5) +
  geom_vline(xintercept=args$min_peak_score, linetype="dashed", color="orange") +
  labs(x="Minimum peak score", y="Total number of peaks") +
  theme(
    axis.text = element_text(size=rel(0.7))
  )

pdf(file.path(args$outdir,"scatterplot_atac_peak_score_vs_number_peaks.pdf"), width = 6, height = 5)
print(p)
dev.off()

#####################################
## Plot fraction of peaks per type ##
#####################################

to.plot <- peak_metadata.dt %>%
  .[score>=args$min_peak_score] %>%
  .[,.N,by="peakType"] %>% .[,percentage:=100*N/sum(N)]

p <- ggpie(to.plot, x="N", label="peakType", fill="peakType") +
  theme(
    legend.position = "none"
  )

pdf(file.path(args$outdir,"pieplot_atac_peak_type.pdf"), width = 5, height = 4)
print(p)
dev.off()

#############################################
## Plot average accessibility per peakType ##
#############################################

to.plot <- peak_metadata.dt %>% 
  .[score>=args$min_peak_score] %>%
  merge(peak_stats.dt,by="peak")

p <- ggboxplot(to.plot, x="peakType", y="mean_cells", fill="peakType", outlier.shape = NA) +
  coord_cartesian(ylim=c(0,1.0)) +
  labs(x="", y="Average chromatin accessibility") +
  theme(
    axis.text = element_text(size=rel(0.75)),
    legend.position = "none"
  )

pdf(file.path(args$outdir,"boxplots_atac_peak_type.pdf"), width = 7, height = 5)
print(p)
dev.off()

# Completion token
file.create(file.path(args$outdir,"completed.txt"))