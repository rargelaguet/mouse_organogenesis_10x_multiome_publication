
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(GGally))
suppressMessages(library(igraph))
suppressMessages(library(network))
suppressMessages(library(sna))
# suppressMessages(library(intergraph))
suppressMessages(library(ggraph))
suppressMessages(library(igraph))
suppressMessages(library(tidygraph))

#####################
## Define settings ##
#####################

# args$rna_cells.sce <- io$rna.sce
io$rna_metacells.sce <- file.path(io$basedir, 'results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds')
io$rna_pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
io$trajectory_file <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacell_trajectory.txt.gz")
io$trajectory <- "nmp"
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacells_metadata.txt.gz")
io$grn_coef <- file.path(io$basedir,'results/rna_atac/gene_regulatory_networks/metacells/trajectories/nmp/global_chip_GRN_coef.txt.gz')
io$outdir <-  file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/metacells/trajectories/nmp/fig"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min_coef <- 0.25
opts$min_tf_score <- 0.75

if (io$trajectory=="nmp") {
  celltypes.to.plot <- c("Caudal_Mesoderm", "Somitic_mesoderm", "NMP", "Spinal_cord")
}

theme_graph <- function() {
  theme_classic() +
    theme(
      axis.line = element_blank(), 
      axis.text = element_blank(),
      axis.ticks = element_blank(), 
      axis.title = element_blank(),
      legend.position = "none"
    )
}

#####################
## Load trajectory ##
#####################

metacell_metadata.dt <- fread(io$metacell_metadata)

trajectory.dt <- fread(io$trajectory_file) %>%
  setnames(c("metacell","FA1","FA2")) %>%
  merge(metacell_metadata.dt[,c("metacell","celltype")])


##############################
## Load marker gene and TFs ##
##############################

# TFs <- fread(paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/TFs.txt"))[[1]]
marker_TFs_all.dt <- fread(io$rna.atlas.marker_TFs.all) %>% 
  .[celltype%in%celltypes.to.plot] %>%
  .[,gene:=toupper(gene)]

marker_TFs_filt.dt <- marker_TFs_all.dt %>% 
  .[score>=opts$min_tf_score]

print(unique(marker_TFs_filt.dt$gene))

marker_TFs_all.dt <- marker_TFs_all.dt[celltype!="Caudal_Mesoderm"]

##############################
## Load RNA expression data ##
##############################

# SingleCellExperiment at cellular resolution
# rna_cells.sce <- readRDS(args$rna_cells.sce)[,trajectory.dt$cell]

# SingleCellExperiment at metacell resolution
rna_metacells.sce <- readRDS(io$rna_metacells.sce)

# SingleCellExperiment at pseudobulk resolution
rna_pseudobulk.sce <- readRDS(io$rna_pseudobulk.sce)#[,celltypes.to.plot]

##################################
## Load global GRN coefficients ##
##################################

GRN_coef.dt <- fread(io$grn_coef) %>% 
  .[,gene:=toupper(gene)] %>% .[gene%in%unique(marker_TFs_filt.dt$gene) & tf%in%unique(marker_TFs_filt.dt$gene)] %>%
  .[pvalue<0.10 & abs(beta)>=opts$min_coef]

# Consider only positive links
# GRN_coef.dt <- GRN_coef.dt[beta>=0]

# GRN_coef.dt[tf=="SOX2" | gene=="T"] %>% View
# GRN_coef.dt[tf=="SOX2"] %>% View

##########################
## Filter TFs and genes ##
##########################

# Filter TFs that have few connections
# TFs <- intersect(GRN_coef.dt[,.N,by="tf"][N>=3,tf], GRN_coef.dt[,.N,by="gene"][N>=2,gene])
TFs <- union(GRN_coef.dt[,.N,by="tf"][N>=3,tf], GRN_coef.dt[,.N,by="gene"][N>=3,gene])
GRN_coef.dt <- GRN_coef.dt[tf%in%TFs & gene%in%TFs]

# Fetch RNA expression matrices
rna_tf_pseudobulk.mtx <- logcounts(rna_pseudobulk.sce)[str_to_title(TFs),]; rownames(rna_tf_pseudobulk.mtx) <- toupper(rownames(rna_tf_pseudobulk.mtx))
rna_tf_metacells.mtx <- logcounts(rna_metacells.sce)[str_to_title(TFs),]; rownames(rna_tf_metacells.mtx) <- toupper(rownames(rna_tf_metacells.mtx))

# Scale
# rna_tf_pseudobulk_scaled.mtx <- apply(rna_tf_pseudobulk.mtx,1,minmax.normalisation) %>% .[celltypes.to.plot,] %>% t 
rna_tf_pseudobulk_scaled.mtx <- apply(rna_tf_pseudobulk.mtx[,celltypes.to.plot],1,minmax.normalisation) %>% t 
rna_tf_metacells_scaled.mtx <- apply(rna_tf_metacells.mtx,1,minmax.normalisation) %>% t 

####################
## Create network ##
####################

# Create node and edge data.frames
TFs <- unique(c(GRN_coef.dt$tf,GRN_coef.dt$gene))
node_list.dt <- data.table(node_id=1:length(TFs), node_name=TFs)
edge_list.dt <- GRN_coef.dt[,c("tf","gene","beta")] %>% setnames(c("from","to","weight")) %>% .[!from==to]

# Create igraph object
igraph.net <- graph_from_data_frame(d = edge_list.dt)

# Create tbl_graph object for ggraph
igraph.tbl <- as_tbl_graph(igraph.net) %>%
  activate(nodes) %>%
  mutate(tf=names(V(igraph.net))) %>%
  mutate(degree=igraph::degree(igraph.net)) %>%
  mutate(eigen_centrality=eigen_centrality(igraph.net)$vector) %>%
  activate(edges) %>%
  mutate(sign=ifelse(E(igraph.net)$weight>0,"Positive","Negative"))

# igraph::closeness(igraph.net)
# igraph::eigen_centrality(igraph.net)$vector
# igraph::betweenness(igraph.net)

#########################
## Plot global network ##
#########################

# Define node color
tmp <- marker_TFs_all.dt[gene%in%names(V(igraph.net))] %>% .[,.SD[which.max(score)][,"celltype"],by="gene"] %>% setkey(gene)
igraph.tbl <- igraph.tbl %>% activate(nodes) %>% mutate(celltype=tmp[names(V(igraph.net)),celltype])

# Plot
set.seed(42)
p <- ggraph(igraph.tbl, 'stress') +
  geom_edge_link(aes(edge_colour = sign), edge_alpha=0.75, edge_width=0.15, arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
  geom_node_point(aes(fill=celltype), size=9, stroke=0, shape=21, alpha=0.75) +
  geom_node_text(aes(label=name), size=3.75) +
  scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_graph()

pdf(file.path(io$outdir,"global_network_celltype.pdf"), width = 5.5, height = 5.75)
print(p)
dev.off()

#################################################
## Plot global network, coloured by centrality ##
#################################################

set.seed(42)
p <- ggraph(igraph.tbl, 'stress') +
  geom_edge_link(edge_colour = "grey66", edge_alpha=0.50, edge_width=0.10, arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
  # geom_node_point(aes(fill=degree), size=8, shape=21) +
  geom_node_point(aes(fill=eigen_centrality), size=8, shape=21, alpha=0.80) +
  geom_node_text(aes(label=name), size=3.75) +
  scale_fill_gradient(low = "white", high = "orange") +
  theme_graph()

pdf(file.path(io$outdir,"global_network_eigen_centrality.pdf"), width = 4.5, height = 5)
print(p)
dev.off()

to.plot <- data.table(tf=names(V(igraph.net)), centrality = eigen_centrality(igraph.net)$vector) %>% 
  .[sample(x=.N, size=round(.N/1.5))] %>%
  setorder(-centrality) %>%
  .[,tf:=factor(tf, levels=tf)] 

# to.plot <- rbind(
#   eigen_centrality.dt %>% head(n=10),
#   eigen_centrality.dt %>% tail(n=10)
# )

p <- ggplot(to.plot, aes_string(x="centrality", y="tf")) +
  geom_point(aes(size=centrality)) +
  geom_segment(aes(yend=tf), xend=0, size=0.20) +
  scale_size_continuous(range = c(1.25,3.75)) +
  labs(y="", x="Eigenvalue centrality") +
  theme_classic() +
  guides(size="none") +
  theme(
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(0.75), color="black")
  )

pdf(file.path(io$outdir,"global_network_eigen_centrality.pdf"), width = 4, height = 4.5)
print(p)
dev.off()

##################################
## Plot repressive interactions ##
##################################

set.seed(42)
p <- ggraph(igraph.tbl, 'stress') +
  geom_edge_link(aes(edge_colour = sign, edge_alpha = sign), edge_width=0.40, arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
  geom_node_point(aes(fill=celltype), size=8, shape=21) +
  geom_node_text(aes(label=name), size=3) +
  scale_edge_alpha_manual(values=c("Positive"=0, "Negative"=0.80)) +
  scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_graph()

pdf(file.path(io$outdir,"global_network_repressive_interactions.pdf"), width = 5, height = 5.5)
print(p)
dev.off()

##################################
## Plot activatory interactions ##
##################################

set.seed(42)
p <- ggraph(igraph.tbl, 'stress') +
  geom_edge_link(aes(edge_colour = sign, edge_alpha = sign), edge_width=0.20, arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
  geom_node_point(aes(fill=celltype), size=8, shape=21) +
  geom_node_text(aes(label=name), size=3) +
  scale_edge_alpha_manual(values=c("Positive"=0.80, "Negative"=0)) +
  scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_graph()

pdf(file.path(io$outdir,"global_network_activatory_interactions.pdf"), width = 5.5, height = 5.5)
print(p)
dev.off()

######################################################################################
## Highlight a TF with pleiotropy effects (both positive and negative interactions) ##
######################################################################################

TF_of_interest <- "CDX2"

idx <- which(names(V(igraph.net))==TF_of_interest)
test <- igraph.tbl %>% activate(edges) %>% filter(from==idx)

tmp <- names(V(igraph.net))[test %>% activate(edges) %>% as_tibble() %>% .$to]
test <- test %>% activate(nodes) %>% filter(name%in%c(TF_of_interest,tmp))

set.seed(42)
p <- ggraph(test, 'stress') +
  geom_edge_link(aes(edge_colour = sign), edge_alpha=0.85, edge_width=0.65, arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
  geom_node_point(aes(fill=celltype), size=11, shape=21, stroke=0, alpha=0.75) +
  geom_node_text(aes(label=name), size=3.5) +
  scale_edge_colour_manual(values=c("Positive"="darkred", "Negative"="darkblue")) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_graph()

pdf(file.path(io$outdir,"CDX2_network.pdf"), width = 4, height = 4)
print(p)
dev.off()

################################
## Plot network per cell type ##
################################

celltypes.to.plot <- c("Spinal_cord","NMP","Somitic_mesoderm")

for (i in celltypes.to.plot) {
  
  # Define node color based on TF expr
  expr.values <- rna_tf_pseudobulk_scaled.mtx[,i]
  expr.values[expr.values<=0.1] <- 0.1; expr.values[expr.values>=0.9] <- 0.9
  igraph.tbl <- igraph.tbl %>% activate(nodes) %>% mutate(expr=expr.values[names(V(igraph.net))])
  
  # (TO-DO) Gray out edges of genes that are not expressed
  
  set.seed(42)
  p <- ggraph(igraph.tbl, 'stress') +
    geom_edge_link(edge_colour="gray70", edge_alpha=0.50, edge_width=0.10, arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
    geom_node_point(aes(color=expr), size=8, stroke=0, alpha=0.75) +
    geom_node_text(aes(label=name), size=3.75) +
    # scale_color_gradient(low = "gray80", high = "#008B00") +
    viridis::scale_color_viridis(option="viridis") +
    theme_graph() + theme(
      legend.position = "right"
    )
  
  pdf(file.path(io$outdir,sprintf("network_coloured_by_%s_expr.pdf",i)), width = 5, height = 5.5)
  print(p)
  dev.off()
  
}

##########
## Plot GIF of the trajectory coloured by expression ##
##########

#######################################
## Plot network along the trajectory ##
#######################################

dir.create(file.path(io$outdir,"gif"), showWarnings = F)

ntimepoints <- 20

trajectory.dt[,pseudotime_group:=as.numeric(cut(FA1,ntimepoints))]

# rna_scaled.mtx <- apply(logcounts(rna_tf_filt.sce),1,minmax.normalisation) %>% t

# i <- "1"
for (i in unique(trajectory.dt$pseudotime_group)) {
  
  cells <- trajectory.dt[pseudotime_group==i,metacell]
  
  rna.expr.i <- apply(rna_tf_metacells_scaled.mtx[,cells],1,mean)
  
  #####################
  ## Plot trajectory ##
  #####################
  
  to.plot <- trajectory.dt %>% copy %>%
    .[,alpha:=1.0] %>% .[!metacell%in%cells,c("celltype","alpha"):=list("None",0.25)]
    
  p1 <- ggplot(to.plot, aes(x=FA1, y=FA2)) +
    geom_point(aes(x=FA1, y=FA2, fill=celltype), size=2.5, shape=21, alpha=1.0, data=to.plot[metacell%in%cells]) +
    geom_point(aes(x=FA1, y=FA2), size=2.5, color="grey", alpha=0.25, data=to.plot[!metacell%in%cells]) +
    # viridis::scale_fill_viridis() +
    scale_fill_manual(values=opts$celltype.colors) +
    # labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position="none"
    )
  
  ##################
  ## Plot network ##
  ##################
  
  # Define node color based on TF expr
  expr.values <- apply(rna_tf_metacells_scaled.mtx[,cells],1,mean)
  expr.values[expr.values<=0.1] <- 0.1; expr.values[expr.values>=0.9] <- 0.9
  igraph.tbl <- igraph.tbl %>% activate(nodes) %>% mutate(expr=expr.values[names(V(igraph.net))])
  
  # (TO-DO) Gray out edges of genes that are not expressed
  tmp <- tibble(tf=names(expr.values), alpha=expr.values, from=1:length(expr.values))
  igraph.tbl.i <- igraph.tbl %>% activate(edges) %>% left_join(tmp,by="from")
  
  
  set.seed(42)
  p2 <- ggraph(igraph.tbl.i, 'stress') +
    # geom_edge_link(edge_colour="gray70", edge_alpha=0.50, edge_width=0.10, arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm')) +
    geom_edge_link(aes(edge_alpha=alpha, edge_width=alpha), edge_colour="gray70", arrow=arrow(length=unit(1.5,'mm')), end_cap=circle(4,'mm'), show.legend = F) +
    scale_edge_alpha_continuous(range=c(0,1)) +
    scale_edge_width_continuous(range=c(0,0.65)) +
    geom_node_point(aes(color=expr), size=14, stroke=0, alpha=0.75) +
    geom_node_text(aes(label=name), size=3.75) +
    # scale_color_gradient(low = "gray80", high = "#008B00") +
    viridis::scale_color_viridis(option="viridis") +
    # guides(edge_alpha="None") +
    theme_graph() + 
    theme(
      legend.position = "right"
    )
  
  
  # Combine plots
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/3,2/3))
  outfile <- sprintf("%s/gif/%d_nmp_pseudotime_GRN.png",io$outdir,i)
  png(outfile, width = 1300, height = 650, bg = "white")
  print(p)
  dev.off()
}


# Create GIF 
library(magick)
sprintf("%s/gif/%d_nmp_pseudotime_GRN.png",io$outdir,sort(unique(trajectory.dt$pseudotime_group))) %>%
  map(image_read) %>%
  image_join %>%
  image_animate(fps=2) %>%
  image_write(quality=100, path = sprintf("%s/gif/test.gif",io$outdir))

