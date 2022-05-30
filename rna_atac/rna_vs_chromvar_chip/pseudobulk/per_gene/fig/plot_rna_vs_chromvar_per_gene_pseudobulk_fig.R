source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

# I/O
io$rna_sce_pseudobulk_file <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_TFs_pseudobulk.rds")
io$atac_chromvar_chip_pseudobulk_file <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_archr.rds",opts$motif_annotation))
io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_gene/fig"); dir.create(io$outdir, showWarnings=F, recursive = T)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load pseudobulk RNA expression
rna_pseudobulk_tf.se <- readRDS(io$rna_sce_pseudobulk_file)

# Load chromVAR matrix
atac_chromvar_pseudobulk.se <- readRDS(io$atac_chromvar_chip_pseudobulk_file)

# Select TFs
TFs <- intersect(rownames(rna_pseudobulk_tf.se),rownames(atac_chromvar_pseudobulk.se))
rna_pseudobulk_tf.se <- rna_pseudobulk_tf.se[TFs,]
atac_chromvar_pseudobulk.se <- atac_chromvar_pseudobulk.se[TFs,]

########################
## Prepare data table ##
########################

atac_chromvar_pseudobulk.dt <- assay(atac_chromvar_pseudobulk.se) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="gene", value.name="chromvar_zscore")

rna_tf_pseudobulk.dt <- logcounts(rna_pseudobulk_tf.se) %>%
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

######################################
## Scatter plot of individual genes ##
######################################

genes.to.plot <- unique(rna_chromvar.dt$gene)
genes.to.plot <- c("FOXA2","FOXB1","FOXC2")

# i <- "FOXA2"
for (i in genes.to.plot) {

  to.plot <- rna_chromvar.dt[gene==i]

  to.plot.text <- rbind(
    to.plot %>% setorder(-expr) %>% head(n=7),
    to.plot %>% setorder(-chromvar_zscore) %>% head(n=7)
  ) %>% unique
  
  p <- ggscatter(to.plot, x="expr", y="chromvar_zscore", fill="celltype", size=5.5, shape=21, 
                  add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
    stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
    ggrepel::geom_text_repel(data=to.plot.text, aes(label=gsub("_"," ",celltype)), size=3.5) +
    scale_fill_manual(values=opts$celltype.colors) +
    # labs(x=sprintf("%s expression",i), y=sprintf("Accessibility of %s targets (z-score)",i)) +
    labs(x="RNA expression", y="chromVAR-ChIP") +
    guides(fill="none") +
    theme(
      axis.text = element_text(size=rel(0.85))
    )
  
  pdf(file.path(io$outdir,sprintf("%s_%s_rna_vs_chromvar_chip_pseudobulk.pdf",i,opts$motif_annotation)), width = 5.5, height = 4)
  print(p)
  dev.off()
}

