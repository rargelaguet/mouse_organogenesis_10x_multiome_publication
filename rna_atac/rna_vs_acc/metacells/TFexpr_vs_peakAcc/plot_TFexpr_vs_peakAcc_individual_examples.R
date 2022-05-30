here::i_am("rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/plot_TFexpr_vs_peakAcc_individual_examples.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='') 
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
p$add_argument('--tf2peak_cor',  type="character",              help='Correlations between TF RNA expression and peak accessibility') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
io$basedir <- file.path(io$basedir,"test")
args <- list()
args$metadata <- file.path(io$basedir,"results/atac/archR/metacells/trajectories/nmp/PeakMatrix/metacells_metadata.txt.gz")
args$sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")
args$atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/metacells/trajectories/nmp/PeakMatrix/PeakMatrix_summarized_experiment_metacells.rds")
args$tf2peak_cor <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/trajectories/nmp/TFexpr_vs_peakAcc/CISBP_cor_TFexpr_vs_peakAcc.rds")
args$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/trajectories/nmp/individual_examples")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)

###################
## Load metadata ##
###################

###############
## Load data ##
###############

# Load SingleCellExperiment
rna_metacells.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_peakMatrix_metacells.se <- readRDS(args$atac_peak_matrix)

# Normalise ATAC data
assayNames(atac_peakMatrix_metacells.se) <- "counts"
assay(atac_peakMatrix_metacells.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_metacells.se),2,colSums(assay(atac_peakMatrix_metacells.se),na.rm=T),"/"))+1)

# Make sure that samples are consistent
samples <- intersect(colnames(rna_metacells.sce),colnames(atac_peakMatrix_metacells.se))
rna_metacells.sce <- rna_metacells.sce[,samples]
atac_peakMatrix_metacells.se <- atac_peakMatrix_metacells.se[,samples]

###############################
## Load TF2peak correlations ##
###############################

tf2peak_cor.se <- readRDS(args$tf2peak_cor)

############
## Filter ##
############

# TFs <- colnames(tf2peak_cor.se)
# peaks <- rownames(tf2peak_cor.se)

##########
## Plot ##
##########

i <- "T"
j <- "chr5:4894870-4895470"

assay(tf2peak_cor.se[j,i],"cor")

stopifnot(colnames(atac_peakMatrix_metacells.se)==colnames(rna_metacells.sce))

to.plot <- data.table(
  atac = assay(atac_peakMatrix_metacells.se[j,],"logcounts")[1,],
  rna = logcounts(rna_metacells.sce[i,])[1,],
  celltype = rna_metacells.sce$celltype
)

p <- ggplot(to.plot, aes(x=rna, y=atac)) +
  geom_point(aes(fill=celltype), color="black", size=2, shape=21) +
  # geom_smooth(method="lm") +
  stat_cor(method = "pearson") +
  scale_fill_manual(values=opts$celltype.colors) +
  # ggrepel::geom_text_repel(aes(label=celltype), size=3, data=to.plot[rna>5 & atac>0.3]) +
  labs(x="RNA expression", y="Peak accessibility", title=sprintf("%s expression vs %s accessibility",i,j)) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust=0.5, size=rel(0.8)),
    axis.text = element_text(color="black"),
    legend.position = "none"
  )

# pdf(file.path(args$outdir,sprintf("%s_vs_%s_rna_vs_acc_metacells.pdf",i,gsub("[:_]","-",j))), width = 8, height = 5)
print(p)
# dev.off()
