here::i_am("rna_atac/rna_vs_chromvar/pseudobulk/per_gene/PAGA/plot_rna_vs_chromvar_paga.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

opts$motif_annotation <- "CISBP"

io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/%s/pdf/paga",opts$motif_annotation))

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
  "Caudal_neurectoderm",
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
genes.to.keep.chromvar <- atac_chromvar_pseudobulk.dt[,.(var(zscore)),by="gene"] %>% .[V1>0.1,gene] %>% as.character
genes.to.keep <- intersect(genes.to.keep.rna, genes.to.keep.chromvar)
rna_tf_pseudobulk.dt <- rna_tf_pseudobulk.dt[gene%in%genes.to.keep]
atac_chromvar_pseudobulk.dt <- atac_chromvar_pseudobulk.dt[gene%in%genes.to.keep]

if (opts$motif_annotation=="CISBP") {
  atac_chromvar_pseudobulk.dt <- atac_chromvar_pseudobulk.dt[!motif=="T_789"]
}

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_pseudobulk.dt,
  by = c("celltype","gene")
)

# Scale values for viz
rna_chromvar.dt[,expr:=minmax.normalisation(expr), by="gene"]
rna_chromvar.dt[,zscore:=minmax.normalisation(zscore), by="gene"]

#####################
## Load PAGA graph ##
#####################

opts$celltypes <- unique(rna_chromvar.dt$celltype)

source(here::here("load_paga_graph.R"))

# Define cell type order
cellype.order <- rownames(connectivity.mtx)

# Parse PAGA graph
igraph.paga.tbl <- igraph.paga.tbl %>% activate(nodes) %>% mutate(diff=diff.values[cellype.order])

#####################################################################
## Plot network, colour by gene expression and motif accessibility ##
#####################################################################

# Define genes to plot
genes.to.plot <- unique(rna_chromvar.dt$gene)# %>% head(n=10)
# genes.to.plot <- c("TAL1","T")

for (i in genes.to.plot) {
  
  outfile <- file.path(io$outdir,sprintf("%s_paga_rna_vs_chromvar.pdf",i))
  if (file.exists(outfile)) {
    print(sprintf("file or %s already exists...",i))
  } else {
    
    # Plot expression
    expr.values <- rna_chromvar.dt[gene==i,expr]; names(expr.values) <- rna_chromvar.dt[gene==i,celltype]
    igraph.paga.tbl <- igraph.paga.tbl %>% activate(nodes) %>% mutate(expr=expr.values[cellype.order])
    
    p1 <- ggraph(igraph.paga.tbl, x = x, y = y) +
      geom_edge_link(edge_colour = "grey66", edge_alpha=1, edge_width=0.5) +
      geom_node_point(aes(fill=expr), size=8, shape=21) +
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
      geom_node_point(aes(fill=chromvar), size=8, shape=21) +
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
    
    pdf(outfile, width = 6.5, height = 4)
    print(p.all)
    dev.off()
  }
}

################################################################################
## (OLD CODE) Plot network, colour by gene expression and motif accessibility ##
################################################################################

# opts$scale <- TRUE
# 
# if (opts$scale) {
#   rna.col.seq <- chromvar.col.seq <- round(seq(0,1,0.1), 2)
#   rna_chromvar.dt[,expr:=minmax.normalisation(expr), by="gene"]
#   rna_chromvar.dt[,zscore:=minmax.normalisation(zscore), by="gene"]
# } else {
#   opts$max.expr <- 9; opts$min.expr <- 2
#   rna.col.seq <- round(seq(opts$min.expr,opts$max.expr,0.1), 2)
#   rna_chromvar.dt[expr>=opts$max.expr,expr:=opts$max.expr]
#   rna_chromvar.dt[expr<=opts$min.expr,expr:=opts$min.expr]
#   
#   opts$max.chromvar <- 6; opts$min.chromvar <- 0
#   chromvar.col.seq <- round(seq(opts$min.chromvar,opts$max.chromvar,0.1), 2)
#   rna_chromvar.dt[zscore>=opts$max.chromvar,zscore:=opts$max.chromvar]
#   rna_chromvar.dt[zscore<=opts$min.chromvar,zscore:=opts$min.chromvar]
#   # activity.col.seq <- round(seq(0,max(rna_chromvar.dt$activity,na.rm=T),0.05), 2)
# }
# 
# # Define colors
# rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))
# chromvar.colors <- colorRampPalette(c("gray92", "purple"))(length(chromvar.col.seq)) 
# # activity.colors <- colorRampPalette(c("gray92", "#EE9A00"))(length(activity.col.seq)) 
# 
# # Define genes to plot
# genes.to.plot <- unique(rna_chromvar.dt$gene)# %>% head(n=10)
# # genes.to.plot <- c("TAL1","T","GLIS2","STAT3")
# 
# # genes.to.plot <- "Prdm1"
# 
# p <- ggnet2(
#   net = net.paga,
#   mode = c("x", "y"),
#   node.size = 0,
#   edge.size = 0.15,
#   edge.color = "grey",
#   label = FALSE,
#   label.size = 2.3
# )
# 
# for (i in genes.to.plot) {
#   
#   if (opts$scale) {
#     outfile <- file.path(io$outdir,sprintf("%s_paga_rna_vs_chromvar_scaled.png",i))
#   } else {
#     outfile <- file.path(io$outdir,sprintf("%s_paga_rna_vs_chromvar.png",i))
#   }
#   
#   if (file.exists(outfile)) {
#     print(sprintf("file or %s already exists...",i))
#   } else {
#     
#     # RNA expression
#     expr.values <- rna_chromvar.dt[gene==i,c("celltype","expr")] %>% matrix.please %>% .[opts$celltypes,]
#     expr.colors <- round(expr.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
#     
#     p1 <- p + geom_text(label = "\u25D0", aes(x=x, y=y), color=expr.colors, size=20, family = "Arial Unicode MS",
#                   data = p$data[,c("x","y")] %>% dplyr::mutate(expr=expr.colors)) +
#       scale_colour_manual(values=expr.colors) + 
#       labs(title="RNA expression") +
#       theme(
#         plot.title = element_text(hjust = 0.5)
#       )
#   
#     # motif accessibility
#     acc.values <- rna_chromvar.dt[gene==i,c("celltype","zscore")] %>% matrix.please %>% .[,1]# %>% .[opts$celltypes,]
#     acc.colors <- round(acc.values,1) %>% map(~ chromvar.colors[which(chromvar.col.seq == .)]) %>% unlist
#     
#     p2 <- p + geom_text(label = "\u25D1", aes(x=x, y=y), color=acc.colors, size=20, family = "Arial Unicode MS",
#                   data = p$data[,c("x","y")] %>% dplyr::mutate(acc=acc.colors)) +
#       scale_fill_manual(values=acc.colors) + 
#       labs(title="Motif accessibility") +
#       theme(
#         plot.title = element_text(hjust = 0.5)
#       )
#     
#     p.all <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1, scale = 0.95) + labs(title=i) + 
#       theme(
#         plot.title = element_text(hjust = 0.5, size=rel(1.25))
#       )
#     
#     png(outfile, width = 650, height = 400)
#     print(p.all)
#     dev.off()
#   }
# }
  