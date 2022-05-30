here::i_am("rna_atac/rna_vs_chromvar_chip/pseudobulk/per_celltype/rna_vs_chromvar_pseudobulk_per_celltype.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--atac_chromvar_chip_pseudobulk',  type="character",              help='ATAC chromVAR-ChIP matrix pseudobulk file') 
p$add_argument('--rna_sce_pseudobulk',  type="character",              help='RNA expression SingleCellExperiment pseudobulk file') 
p$add_argument('--tf_markers_rna',  type="character",              help='') 
p$add_argument('--tf_markers_atac',  type="character",              help='') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args <- list()
args$motif_annotation <- "CISBP"
args$rna_sce_pseudobulk <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_TFs_pseudobulk.rds")
args$atac_chromvar_chip_pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_archr.rds",args$motif_annotation))
# args$atac_chromvar_chip_pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_chromvar.rds",args$motif_annotation))
args$tf_markers_rna <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_tfs_all.txt.gz")
args$tf_markers_atac <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltypes/CISBP/parsed/markers_all.txt.gz")
args$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_celltype/%s",args$motif_annotation))
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
dir.create(args$outdir, showWarnings=F, recursive = T)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load pseudobulk TF RNA expression
rna_pseudobulk.sce <- readRDS(args$rna_sce_pseudobulk)

# Load pseudobulk ATAC chromVAR-ChIP matrix
atac_pseudobulk_chromvar.se <- readRDS(args$atac_chromvar_chip_pseudobulk)

# Fetch common TFs
TFs <- intersect(rownames(rna_pseudobulk.sce),rownames(atac_pseudobulk_chromvar.se))
rna_pseudobulk.sce <- rna_pseudobulk.sce[TFs,]
atac_pseudobulk_chromvar.se <- atac_pseudobulk_chromvar.se[TFs,]

########################
## Prepare data table ##
########################

atac_chromvar_pseudobulk.dt <- assay(atac_pseudobulk_chromvar.se,"z") %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="gene", value.name="chromvar_zscore")

rna_tf_pseudobulk.dt <- logcounts(rna_pseudobulk.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  data.table::melt(id.vars="gene", variable.name="celltype", value.name="expr")

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_pseudobulk.dt,
  by = c("celltype","gene")
)

length(unique(rna_chromvar.dt$gene))

#############################################
## Load pre-computed correlation estimates ##
#############################################

# cor.dt <- fread(args$cor_rna_vs_chromvar_per_gene) %>%
#   .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])]

###########################
## Load TF marker scores ##
###########################

tf_markers_rna.dt <- fread(args$tf_markers_rna) %>% setnames("score","tf_marker_score_rna")
tf_markers_atac.dt <- fread(args$tf_markers_atac) %>% setnames("score","tf_marker_score_atac")

###############################
## Scatterplot per cell type ##
###############################

# opts$max.chromvar <- 12
# opts$max.expr <- 12
# opts$min.expr <- 3

celltypes.to.plot <- unique(rna_chromvar.dt$celltype)

# i <- "Neural_crest"
for (i in celltypes.to.plot) {
  
  to.plot <- rna_chromvar.dt[celltype==i]  %>%
    merge(tf_markers_rna.dt[celltype==i,c("celltype","gene","tf_marker_score_rna")],by=c("gene","celltype")) %>%
    merge(tf_markers_atac.dt[celltype==i,c("celltype","gene","tf_marker_score_atac")],by=c("gene","celltype")) %>%
    # merge(cor.dt[,c("gene","cor_sign")], by="gene") %>%
    .[chromvar_zscore<0,chromvar_zscore:=0] %>%
    # .[chromvar_zscore>=opts$max.chromvar,chromvar_zscore:=opts$max.chromvar] %>%
    .[expr>=opts$max.expr,expr:=opts$max.expr] %>%
    .[expr<=opts$min.expr,expr:=opts$min.expr]
  
  # Remove TFs that are not expressed in this cell type
  to.plot <- to.plot[expr>=2.5]
  
  to.plot %>% .[,dot_size:=minmax.normalisation(tf_marker_score_rna*tf_marker_score_atac)]
  to.plot.text <- to.plot[tf_marker_score_rna>=0.85 & tf_marker_score_atac>=0.85]
  # to.plot.dots <- to.plot[!gene%in%to.plot.text$gene]
  
  
  p <- ggplot(to.plot, aes(x=tf_marker_score_atac, y=tf_marker_score_rna)) +
    geom_jitter(aes(size=dot_size), shape=21, width=0.01, height=0.01, fill="gray80") + 
    # geom_point(size=0.5, data=to.plot.text) +
    ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=3, max.overlaps=Inf) +
    # geom_text(aes(label=gene), size=3, data=to.plot.text) +
    scale_size_continuous(range = c(0.1,4)) +
    scale_fill_gradient(low = "gray80", high = "darkgreen") +
    # scale_fill_brewer(palette="Dark2") +
    # coord_cartesian(
    #   # xlim = c(0,opts$max.chromvar+0.1),
    #   ylim = c(opts$min.expr,opts$max.expr+0.20)
    #   ) +
    labs(x="TF marker score (ATAC)", y="TF marker score (RNA)") +
    # guides(size="none", fill = guide_legend(override.aes = list(size=3))) +
    guides(size="none") +
    theme_classic() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.text = element_text(size=rel(0.75), color="black")
    )
  
  outfile <- file.path(args$outdir,sprintf("%s_rna_vs_chromvar_chip_%s_pseudobulk.pdf",i,args$motif_annotation))
  pdf(outfile, width = 8, height = 6)
  print(p)
  dev.off()
} 

# Create a completion token
file.create(file.path(args$outdir,sprintf("%s_completed.txt",args$motif_annotation)))

###############################
## Polar plot per cell type ##
###############################

celltypes.to.plot <- unique(rna_chromvar.dt$celltype)

# i <- "Neural_crest"
for (i in celltypes.to.plot) {
  
  to.plot <- rna_chromvar.dt[celltype==i]  %>%
    merge(tf_markers_rna.dt[celltype==i,c("celltype","gene","tf_marker_score_rna")],by=c("gene","celltype")) %>%
    merge(tf_markers_atac.dt[celltype==i,c("celltype","gene","tf_marker_score_atac")],by=c("gene","celltype"))
  
  to.plot %>% .[,tf_marker_score:=tf_marker_score_rna*tf_marker_score_rna] %>% .[,tf_marker_score:=tf_marker_score/max(tf_marker_score)]
  to.plot %>% .[,dot_size:=minmax.normalisation(tf_marker_score)]
  
  to.plot.text <- to.plot[tf_marker_score_rna>=0.85 & tf_marker_score_atac>=0.85]
  
  p <- ggplot(to.plot, aes(x=gene, y=tf_marker_score, group=1)) +
    geom_point(aes(size=dot_size), stat = 'identity', shape=21, color="black", fill="gray95", alpha=0.9) +
    ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=3, max.overlaps=Inf) +
    # geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line[celltype%in%to.plot$celltype]) +
    scale_size_continuous(range = c(0.05,3)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by=0.15)) +
    guides(size="none") +
    coord_polar() + 
    theme_bw() +
    theme(
      panel.background=element_rect(fill = c("white")),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(size=rel(1.5)),
      axis.title=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_text(),
      axis.ticks.y=element_blank(),
      axis.line.x=element_blank()
    )
  
  outfile <- file.path(args$outdir,sprintf("%s_rna_vs_chromvar_chip_%s_polar_plots_pseudobulk.pdf",i,args$motif_annotation))
  pdf(outfile, width = 8, height = 6)
  print(p)
  dev.off()
  
}

##########
## TEST ##
##########

# to.plot <- rna_chromvar.dt[gene%in%c("SOX2","T")] %>% dcast(celltype~gene, value.var=c("expr","chromvar_zscore"))
# 
# ggscatter(to.plot, x="chromvar_zscore_SOX2", y="chromvar_zscore_T", fill="celltype", shape=21, size=4) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   theme(
#     legend.position = "none"
#   )
