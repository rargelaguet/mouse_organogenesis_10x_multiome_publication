here::i_am("rna_atac/rna_vs_chromvar/pseudobulk/per_gene/PAGA/plot_rna_vs_chromvar_paga.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

opts$motif_annotation <- "CISBP"

io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_gene/%s/pdf/paga/fig",opts$motif_annotation))
dir.create(io$outdir, showWarnings = F, recursive = T)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

################################################
## Load pseudobulk RNA and chromVAR estimates ##
################################################

source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

################
## Parse data ##
################

# Filter genes with low variability
genes.to.keep.rna <- rna_tf_pseudobulk.dt[,.(var(expr)),by="gene"] %>% .[V1>0.1,gene] %>% as.character
genes.to.keep.chromvar <- atac_chromvar_chip_pseudobulk.dt[,.(var(zscore)),by="gene"] %>% .[V1>0.1,gene] %>% as.character
genes.to.keep <- intersect(genes.to.keep.rna, genes.to.keep.chromvar)
rna_tf_pseudobulk.dt <- rna_tf_pseudobulk.dt[gene%in%genes.to.keep]
atac_chromvar_chip_pseudobulk.dt <- atac_chromvar_chip_pseudobulk.dt[gene%in%genes.to.keep]

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_chip_pseudobulk.dt,
  by = c("celltype","gene")
)

#####################
## Load PAGA graph ##
#####################

opts$celltypes <- unique(rna_chromvar.dt$celltype)

source(here::here("load_paga_graph.R"))

# Define cell type order
cellype.order <- rownames(connectivity.mtx)

#####################################################################
## Plot network, colour by gene expression and motif accessibility ##
#####################################################################

# Define genes to plot
genes.to.plot <- c("FOXA2","TAL1","RFX4")
genes.to.plot <- c("ETS1","TFAP2C","SOX9")

for (i in genes.to.plot) {

  # Plot expression
  expr.values <- rna_chromvar.dt[gene==i,expr]; names(expr.values) <- rna_chromvar.dt[gene==i,celltype]
  expr.values[expr.values<=1.5] <- 1.5
  
  igraph.paga.tbl <- igraph.paga.tbl %>% activate(nodes) %>% mutate(expr=expr.values[cellype.order])
  
  p1 <- ggraph(igraph.paga.tbl, x = x, y = y) +
    geom_edge_link(edge_colour = "grey66", edge_alpha=1, edge_width=0.5) +
    geom_node_point(aes(fill=expr), size=6, shape=21) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_classic(base_size=14) +
    theme(
      axis.line = element_blank(), 
      axis.text = element_blank(),
      axis.ticks = element_blank(), 
      axis.title = element_blank(),
      legend.position = "none"
    )
  
  # Plot chromvar
  chromvar.values <- rna_chromvar.dt[gene==i,zscore]; names(chromvar.values) <- rna_chromvar.dt[gene==i,celltype]
  igraph.paga.tbl <- igraph.paga.tbl %>% activate(nodes) %>% mutate(chromvar=chromvar.values[cellype.order])
  
  p2 <- ggraph(igraph.paga.tbl, x = x, y = y) +
    geom_edge_link(edge_colour = "grey66", edge_alpha=1, edge_width=0.5) +
    geom_node_point(aes(fill=chromvar), size=6, shape=21) +
    scale_fill_gradient(low = "white", high = "purple") +
    theme_classic(base_size=14) +
    theme(
      axis.line = element_blank(), 
      axis.text = element_blank(),
      axis.ticks = element_blank(), 
      axis.title = element_blank(),
      legend.position = "none"
    )
  
  p.all <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1, scale = 0.95) + labs(title=i) + 
    theme(
      plot.title = element_text(hjust = 0.5, size=rel(1.25))
    )
    
    pdf(file.path(io$outdir,sprintf("%s_paga_rna_vs_chromvar.pdf",i)), width = 6, height = 3.25)
    print(p.all)
    dev.off()
}


