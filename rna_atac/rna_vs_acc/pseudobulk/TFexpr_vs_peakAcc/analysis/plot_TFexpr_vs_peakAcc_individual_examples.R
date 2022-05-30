here::i_am("rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/plot_TFexpr_vs_peakAcc_individual_examples.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--atac_peak_matrix',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
p$add_argument('--tf2peak_cor',  type="character",              help='Correlations between TF RNA expression and peak accessibility') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
# args$atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype/PeakMatrix/pseudobulk_PeakMatrix_summarized_experiment.rds")
# args$tf2peak_cor <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/CISBP_cor_TFexpr_vs_peakAcc.rds")
# args$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/individual_examples")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)

###############
## Load data ##
###############

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(args$sce)

# Load ATAC SummarizedExperiment
atac_peakMatrix_pseudobulk.se <- readRDS(args$atac_peak_matrix)

# Normalise ATAC data
assayNames(atac_peakMatrix_pseudobulk.se) <- "counts"
assay(atac_peakMatrix_pseudobulk.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_pseudobulk.se),2,colSums(assay(atac_peakMatrix_pseudobulk.se),na.rm=T),"/"))+1)

###############################
## Load TF2peak correlations ##
###############################

tf2peak_cor.se <- readRDS(args$tf2peak_cor)

############
## Filter ##
############

TFs <- colnames(tf2peak_cor.se)
peaks <- rownames(tf2peak_cor.se)

##########
## Plot ##
##########

# tmp <- assay(tf2peak_cor.se[,i],"cor")[,1]
i <- "T"
j <- "chr5:4894870-4895470"

to.plot <- data.table(
  atac = assay(atac_peakMatrix_pseudobulk.se[j,],"logcounts")[1,],
  rna = logcounts(rna_pseudobulk.sce[i,])[1,],
  celltype = colnames(rna_pseudobulk.sce)
)

p <- ggplot(to.plot, aes(x=rna, y=atac, fill=celltype)) +
  geom_point(color="black", size=4, shape=21) +
  # geom_smooth(method="lm") +
  stat_cor(method = "pearson") +
  scale_fill_manual(values=opts$celltype.colors) +
  ggrepel::geom_text_repel(aes(label=celltype), size=3, data=to.plot[rna>5 & atac>0.3]) +
  labs(x="RNA expression", y="Peak accessibility", title=sprintf("%s expression vs %s accessibility",i,j)) +
  theme_classic() + 
  theme(
    plot.title = element_text(hjust=0.5, size=rel(0.8)),
    axis.text = element_text(color="black"),
    legend.position = "none"
  )

pdf(file.path(args$outdir,sprintf("%s_vs_%s_rna_vs_acc_pseudobulk.pdf",i,gsub("[:_]","-",j))), width = 8, height = 5)
print(p)
dev.off()