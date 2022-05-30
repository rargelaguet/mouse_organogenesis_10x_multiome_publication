# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$pseudobulk_atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype_genotype/PeakMatrix/PeakMatrix_pseudobulk_with_replicates.rds")
io$outdir <- file.path(io$basedir,"results/atac/archR/plot_individual_peaks/genotype"); dir.create(io$outdir, showWarnings = F)

# Options
opts$samples <- c(
  "E8.5_CRISPR_T_KO",
  "E8.5_CRISPR_T_WT"
)

opts$celltypes <- c("Somitic_mesoderm", "NMP", "Spinal_cord")

####################
## Load metadata  ##
####################

##########################
## Load ATAC PeakMatrix ##
##########################

atac_peak_matrix_pseudobulk.se <- readRDS(io$pseudobulk_atac_peak_matrix)

# subset
atac_peak_matrix_pseudobulk.se <- atac_peak_matrix_pseudobulk.se[,atac_peak_matrix_pseudobulk.se$celltype%in%opts$celltypes]

# Normalise ATAC data
assay(atac_peak_matrix_pseudobulk.se,"logcounts") <- log(1e6*(sweep(assay(atac_peak_matrix_pseudobulk.se),2,colSums(assay(atac_peak_matrix_pseudobulk.se),na.rm=T),"/"))+1)

####################################################
## Boxplots of chromatin accessibility (WT vs KO) ##
####################################################

peaks.to.plot <- c("chr7:144884955-144885555","chr7:79789147-79789747","chr7:126785067-126785667")

# i <- "chr7:144884955-144885555"
for (i in peaks.to.plot) {
  
  to.plot <- data.table(
    acc = assay(atac_peak_matrix_pseudobulk.se,"logcounts")[i,],
    sample = colnames(atac_peak_matrix_pseudobulk.se),
    celltype = atac_peak_matrix_pseudobulk.se$celltype,
    genotype = atac_peak_matrix_pseudobulk.se$genotype
  ) %>% .[celltype=="Caudal_Mesoderm",celltype:="Somitic_mesoderm"] %>% 
    .[,celltype_genotype:=sprintf("%s (%s)",celltype,genotype)]
  
  order <- c("Spinal_cord (WT)","Spinal_cord (T_KO)", "NMP (WT)", "NMP (T_KO)", "Somitic_mesoderm (WT)")
  to.plot[,celltype_genotype:=factor(celltype_genotype, levels=order)]
  
  my_comparisons <- list( c("NMP (WT)", "NMP (T_KO)"))
  
  to.plot.means <- to.plot[,.(acc=mean(acc),sd=sd(acc)), by=c("celltype_genotype","celltype","genotype")]
  
  p <- ggplot(to.plot, aes_string(x="celltype_genotype", y="acc", fill="genotype")) +
    geom_bar(stat="identity", color="black", alpha=1, data=to.plot.means) +
    geom_jitter(size=3, width=0.05, shape=21) +
    geom_errorbar(aes(ymin=acc-sd, ymax=acc+sd), width=0.15, alpha=1, size=0.6, data=to.plot.means) +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), comparisons = my_comparisons, method="t.test") +
    # stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width=0.75)) +
    scale_fill_manual(values=c("#EE0000","#1C86EE")) +
    labs(x="", y="Chromatin accessibility (log normalised counts)") +
    # geom_violin(aes(fill=celltype)) +
    theme_classic() +
    theme(
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black"),
      legend.title = element_blank(),
      legend.position = "none"
    )
  
  pdf(file.path(io$outdir,sprintf("boxplots_acc_genotype_%s.pdf",gsub(":","-",i))), width=6, height=6)
  print(p)
  dev.off()
}
