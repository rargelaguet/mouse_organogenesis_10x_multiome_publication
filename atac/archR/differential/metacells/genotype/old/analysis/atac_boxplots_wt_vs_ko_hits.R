# here::i_am("atac/archR/differential/wt_vs_ko/analysis/analysis.R")

source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("atac/archR/differential/utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$min.MeanDiff <- 0.10
opts$threshold_fdr <- 0.01
opts$atac.matrix <- "PeakMatrix"

opts$wt.class <- "WT"
opts$ko.class <- "KO"

# I/O
io$diff.dir <- file.path(io$basedir,sprintf("results/atac/archR/differential/wt_vs_ko/%s",opts$atac.matrix))
io$atac_pseudobulk_peak_matrix.se <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype_genotype/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/wt_vs_ko/%s/pdf/individual_hits",opts$atac.matrix)); dir.create(io$outdir, showWarnings = F)

##########################
## Load pseudobulk ATAC ##
##########################

atac_pseudobulk_peakMatrix.se <- readRDS(io$atac_pseudobulk_peak_matrix.se)

atac_pseudobulk_peakMatrix.se$celltype <- strsplit(colnames(atac_pseudobulk_peakMatrix.se),split = "-") %>% map_chr(1)
atac_pseudobulk_peakMatrix.se$genotype <- strsplit(colnames(atac_pseudobulk_peakMatrix.se),split = "-") %>% map_chr(2)

celltypes.to.plot <- which(table(atac_pseudobulk_peakMatrix.se$celltype)==2) %>% names

io$atac_pseudobulk_stats <- "/Users/argelagr/data/gastrulation_multiome_10x/results/atac/archR/pseudobulk/celltype_genotype/stats.txt"
atac_pseudobulk_stats.dt <- fread(io$atac_pseudobulk_stats) %>%
  .[,celltype:=strsplit(group,split = "-") %>% map_chr(1)] %>%
  .[,genotype:=strsplit(group,split = "-") %>% map_chr(2)]
celltypes.to.plot <- opts$celltypes[opts$celltypes%in%names(which(table(atac_pseudobulk_stats.dt[N>=75,celltype])==2))]

celltypes.to.plot <- setdiff(celltypes.to.plot,c("Erythroid1","ExE_endoderm","ExE_mesoderm","Parietal_endoderm","Pharyngeal_mesoderm"))

###############################
## Load differential results ##
###############################

diff.dt <- opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_WT_vs_KO.txt.gz", io$diff.dir,j)
  if (file.exists(file)) {
    fread(file, select=c(1,2,3)) %>% .[,c("celltype"):=list(j)]
  }
}) %>% rbindlist %>%
  .[, sign := ifelse(MeanDiff>0,"Downregulated in KO","Upregulated in KO")] %>%
  .[, sig := (FDR<=opts$threshold_fdr & abs(MeanDiff)>=opts$min.MeanDiff)]


##############
## Boxplots ##
##############

peaks.to.plot <- diff.dt[celltype=="NMP" & sign=="Downregulated in KO" & abs(MeanDiff)>=0.15 & sig==TRUE,idx]

# i <- peaks.to.plot[1]
for (i  in peaks.to.plot) {
  
  to.plot <- data.table(
    celltype = atac_pseudobulk_peakMatrix.se$celltype,
    genotype = atac_pseudobulk_peakMatrix.se$genotype,
    acc = assay(atac_pseudobulk_peakMatrix.se[i,])[1,]
  ) %>% .[celltype%in%celltypes.to.plot]
  
  sig.dt <- diff.dt[idx==i & celltype%in%celltypes.to.plot,c("celltype","FDR")] %>% 
    .[,label:=sprintf("p=%s",round(FDR,2))] %>%
    .[,y:=max(to.plot$acc)+0.03] 
  
  p <- ggplot(to.plot, aes(x=celltype, y=acc)) +
    geom_bar(aes(fill = genotype), color="black", stat = 'identity', position="dodge") + 
    geom_text(aes(y=y, label=label), data=sig.dt, size=3) + 
    # scale_fill_manual(values=opts$celltype.colors, drop=F) +
    guides(x = guide_axis(angle = 90)) +
    labs(x="", y="Chromatin accessibility levels", title=i) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size=rel(0.75), color="black"),
      axis.text.y = element_text(size=rel(0.9), color="black")
    )
  
  pdf(file.path(io$outdir,sprintf("%s_wt_vs_ko_barplots_acc.pdf",gsub("[:_]","-",i))), width=9, height=6)
  print(p)
  dev.off()
}
