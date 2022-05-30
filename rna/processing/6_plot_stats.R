here::i_am("rna/processing/6_plot_stats.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--samples',       type="character",  default="all",  nargs='+',  help='Samples to plot')
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args <- list()
args$metadata <- file.path(io$basedir,"results/rna/doublet_detection/sample_metadata_after_doublets.txt.gz")
args$sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
args$samples <- c("E8.5_CRISPR_T_KO", "E8.5_CRISPR_T_WT") # "all"
args$outdir <- file.path(io$basedir,"results/rna/stats")
## END TEST ##

#####################
## Parse arguments ##
#####################

# I/O
dir.create(args$outdir, showWarnings = T)

# Options
if (args$samples[1]=="all") {
  args$samples <- opts$samples
} else {
  stopifnot(args$samples%in%opts$samples)
}

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$samples] %>%
  .[,sample:=factor(sample,levels=opts$samples[opts$samples%in%args$samples])]

# sample_metadata <- fread(io$metadata) %>%  .[(pass_rnaQC==TRUE | pass_atacQC==TRUE)]
# table(sample_metadata$sample)
      
###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

################################################
## Plot distribution of UMI counts per sample ##
################################################

to.plot <- args$samples %>% map(function(i) {
  mtx <- as.matrix(counts(sce[,sce$sample==i])[,1:100])
  data.table(table(mtx[mtx<=10])) %>% setnames(c("value","N")) %>% .[,sample:=i]
}) %>% rbindlist

to.plot[N>=1e6,N:=1e6]

p <- ggbarplot(to.plot, x="value", y="N", fill="gray70") +
  facet_wrap(~sample) +
  theme(
    axis.text =  element_text(size=rel(0.65)),
    axis.title.x = element_blank()
  )

pdf(paste0(args$outdir,"/umi_distribution_per_sample.pdf"), width=8, height=5)
print(p)
dev.off()

#########################################
## Barplots number of cells per sample ##
#########################################

to.plot <- sample_metadata[,.N,by=c("sample","stage")]
fwrite(to.plot[,c("sample","N")], file.path(args$outdir,"ncells_per_sample.txt.gz"), sep="\t", quote=F)

p <- ggbarplot(to.plot, x = "sample", y = "N", fill="stage") +
  scale_fill_manual(values=opts$stage.colors) +
  labs(x="", y="Number of cells (after QC)") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(0.7)),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=40, hjust=1)
  )

pdf(paste0(args$outdir,"/ncells_per_sample.pdf"), width=8, height=5)
print(p)
dev.off()

#########################################
## Barplots number of reads per sample ##
#########################################

to.plot <- data.table(cell=colnames(sce), sample=sce$sample, stage=sce$stage, N=colSums(counts(sce))) %>%
  .[,.(N=sum(N)),by=c("sample","stage")]
fwrite(to.plot[,c("sample","N")], file.path(args$outdir,"nreads_per_sample.txt.gz"), sep="\t", quote=F)

p <- ggbarplot(to.plot, x = "sample", y = "N", fill="stage") +
  scale_fill_manual(values=opts$stage.colors) +
  labs(x="", y="Number of reads (after QC)") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(0.7)),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=40, hjust=1)
  )

pdf(paste0(args$outdir,"/ncells_per_sample.pdf"), width=8, height=5)
print(p)
dev.off()

