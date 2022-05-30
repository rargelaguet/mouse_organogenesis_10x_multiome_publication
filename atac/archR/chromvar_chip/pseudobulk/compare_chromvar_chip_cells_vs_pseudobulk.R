here::i_am("atac/archR/chromvar_chip/compare_chromvar_chip_cells_vs_pseudobulk.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$chromvar_chip_cells <- file.path(io$basedir,"results_new/atac/archR/chromvar_chip/chromVAR_chip_deviations_CISBP_archr.rds")
io$chromvar_chip_pseudobulk <- file.path(io$basedir,"results_new/atac/archR/chromvar_chip/pseudobulk/chromVAR_deviations_CISBP_archr_chip.rds")
io$outdir <- file.path(io$basedir,"results_new/atac/archR/chromvar_chip/pseudobulk/comparison"); dir.create(io$outdir, showWarnings = F)

# Options
opts$motif_annotation <- "CISBP"

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[,celltype:=factor(celltype.predicted,levels=names(opts$celltype.colors))]

################################################################
## Fetch chromVAR-ChIP deviations computed using single cells ##
################################################################

chromvar_deviations_cells.se <- readRDS(io$chromvar_chip_cells)

##############################################################
## Fetch chromVAR-ChIP deviations computed using pseudobulk ##
##############################################################

chromvar_deviations_pseudobulk.se <- readRDS(io$chromvar_chip_pseudobulk)

#############
## Combine ##
#############

# chromvar.deviations.dt <- rbind(
#   chromvar_deviations_cells.dt,
#   chromvar_deviations_pseudobulk.dt
# )

# length(unique(chromvar.deviations.dt$motif))
# length(unique(chromvar.deviations.dt$cell))

##########
## Plot ##
##########

# genes.to.plot <- c("GATA1","FOXA2","TFAP2A")
genes.to.plot <- intersect(rownames(chromvar_deviations_cells.se), rownames(chromvar_deviations_pseudobulk.se))

opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(sample_metadata$celltype)]

opts$ymin_pseudobulk <- -3 # min(assay(chromvar_deviations_pseudobulk.se))
opts$ymax_pseudobulk <- 4 # max(assay(chromvar_deviations_pseudobulk.se))
opts$min_value_pseudobulk <- -3.1
opts$max_value_pseudobulk <- 4.1

# opts$ymin_cells <- min(assay(chromvar_deviations_cells.se),na.rm=T)
# opts$ymax_cells <- max(assay(chromvar_deviations_cells.se),na.rm=T)
opts$ymin_cells <- -3
opts$ymax_cells <- 11
opts$min_value_cells <- -3.1
opts$max_value_cells <- 11.1

# i <- "GATA1"
for (i in genes.to.plot) {
  
  pseudobulk_to_plot <- data.table(
    celltype = colnames(chromvar_deviations_pseudobulk.se),
    value = assay(chromvar_deviations_pseudobulk.se[i,])[1,]
  ) %>% .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))] %>%
    .[value>opts$max_value_pseudobulk,value:=opts$max_value_pseudobulk] %>%
    .[value<opts$min_value_pseudobulk,value:=opts$min_value_pseudobulk]

  p.pseudobulk <- ggbarplot(pseudobulk_to_plot, x="celltype", y="value", fill="celltype") +
    scale_fill_manual(values=opts$celltype.colors) +
    geom_hline(yintercept=0, linetype="dashed") +
    coord_cartesian(ylim=c(opts$ymin_pseudobulk,opts$ymax_pseudobulk)) +
    guides(x = guide_axis(angle = 90)) +
    labs(x="", y=sprintf("%s chromVAR-ChIP z-score",i)) +
    theme(
      legend.position = "none",
      # axis.text.x = element_text(color="black", angle=40, hjust=1, size=rel(0.75)),
      axis.text = element_text(color="black", size=rel(0.6)),
      axis.title = element_text(color="black", size=rel(0.7)),
      axis.ticks.x = element_blank(),
    )

  cells_to_plot <- data.table(
    cell = colnames(chromvar_deviations_cells.se),
    value = assay(chromvar_deviations_cells.se[i,])[1,]
  ) %>% merge(sample_metadata[,c("cell","celltype")], by="cell") %>%
    .[value>opts$max_value_cells,value:=opts$max_value_cells] %>%
    .[value<opts$min_value_cells,value:=opts$min_value_cells]

  p.cells <- ggboxplot(cells_to_plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
    geom_hline(yintercept=0, linetype="dashed") +
    coord_cartesian(ylim=c(opts$ymin_cells,opts$ymax_cells)) +
    guides(x = guide_axis(angle = 90)) +
    labs(x="", y=sprintf("%s chromVAR-ChIP z-score",i)) +
    theme(
      legend.position = "none",
      axis.text = element_text(color="black", size=rel(0.6)),
      axis.title = element_text(color="black", size=rel(0.7)),
      axis.ticks.x = element_blank()
  
    )

  p <- cowplot::plot_grid(plotlist=list(p.pseudobulk,p.cells), nrow=1, rel_widths = c(1/2,1/2))

  pdf(file.path(io$outdir,sprintf("%s_chromvar_chip_cells_vs_pseudobulk.pdf",i)), width=10, height=4)
  print(p)
  dev.off()
}
