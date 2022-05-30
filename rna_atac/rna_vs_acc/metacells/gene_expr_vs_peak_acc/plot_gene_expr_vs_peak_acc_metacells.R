here::i_am("rna_atac/rna_vs_acc/pseudobulk/gene_expr_vs_promoter_acc/plot_gene_expr_vs_promoter_acc_pseudobulk.R")


source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

stop()
p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (pseudobulk)') 
p$add_argument('--gene_score_matrix',  type="character",              help='ATAC Gene score matrix (pseudobulk)') 
p$add_argument('--outfile',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args <- list()
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped.rds") 
# args$gene_score_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype.mapped/pseudobulk_GeneScoreMatrix_TSS_summarized_experiment.rds") # io$archR.pseudobulk.GeneMatrix.se
# args$outfile <- file.path(io$basedir,"results/rna_atac/gene_expr_vs_promoter_acc/pseudobulk/cor_gene_expr_vs_promoter_acc_pseudobulk.txt.gz")
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
io$gene_expr_vs_peak_acc <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/gene_expr_vs_peak_acc/cor_gene_expr_vs_peak_acc_nearest.txt.gz")
io$outdir <- file.path(io$basedir,"results_new/rna_atac/rna_vs_acc/pseudobulk/gene_expr_vs_peak_acc/pdf")

# Options

###############
## Load data ##
###############

gene_expr_vs_peak_acc.dt <- fread(io$gene_expr_vs_peak_acc) %>%
  .[is.na(cor), c("cor","pvalue"):=list(0,1)]

# Filter genes
# grep("Rik",unique(gene_expr_vs_peak_acc.dt$gene), value=T)

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(io$archR.peak.metadata) %>%
  # .[peakType!="Promoter"] %>%
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)] 


# For promoter peaks, select the one with the shortest distance to the TSS
peak_metadata.dt <- rbind(
  peak_metadata.dt[peakType=="Promoter" & distToTSS<=200,.SD[which.min(distToTSS)],by="nearestGene"],
  peak_metadata.dt[peakType!="Promoter"]
)

stopifnot(peak_metadata.dt[,.N,by=c("peak","nearestGene")][,N]==1)
# stopifnot(unique(gene_expr_vs_peak_acc.dt$peak)%in%peak_metadata.dt$peak)

gene_expr_vs_peak_acc.dt <- gene_expr_vs_peak_acc.dt %>% 
  merge(peak_metadata.dt[,c("peak","peakType","distToTSS","nearestGene")],by="peak") %>%
  .[,sign:=c("-","+")[as.numeric(cor>0)+1]]

###########################################################################################
## Calculate number of peaks that are correlated with RNA expression of the nearest gene ##
###########################################################################################

to.plot <- gene_expr_vs_peak_acc.dt %>%
  .[,sum(pvalue<0.05 & abs(cor)>=0.25), by=c("peak","sign","peakType")] %>%
  .[,mean(V1>=1),by=c("sign","peakType")]

# ggbarplot(to.plot, x="peakType", fill="sign", y="V1", position=position_dodge(width = 0.75)) +
ggbarplot(to.plot, x="peakType", fill="sign", y="V1") +
  labs(x="", y="Fraction of peaks correlated with gene expr.") +
  # scale_fill_brewer(palette="Dark2") +
  # guides(color=F) + 
  theme(
    axis.title = element_text(size=rel(0.80)),
    legend.title = element_blank(),
    axis.text.x = element_text(size=rel(0.80), color="black"),
    axis.text.y = element_text(size=rel(0.75), color="black")
  )

##################################################
## Boxplot of ATAC peak acc vs RNA correlations ##
##################################################

to.plot <- gene_expr_vs_peak_acc.dt %>% .[pvalue<=0.05]

ggboxplot(to.plot[sign=="+"], x="peakType", y="cor") +
  # facet_wrap(~sign) +
  # geom_vline(aes(xintercept=V1), color="black", data=to.plot[,median(value),by="variable"], linetype="dashed") +
  # labs(x="Number of TF motifs per peak", y="Density") +
  # scale_x_continuous(limits = c(0,165)) +
  scale_fill_brewer(palette="Dark2") +
  # guides(color=F) +
  theme(
    axis.title = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.text.x = element_text(size=rel(0.80), color="black"),
    axis.text.y = element_text(size=rel(0.50), color="black")
  )


pdf(file.path(io$outdir,"density_number_motifs_per_peak_before_and_after.pdf"), width=5, height=4)
print(p)
dev.off()
  