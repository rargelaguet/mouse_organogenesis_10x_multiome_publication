# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

## I/O
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacells_metadata.txt.gz")
io$metacell_sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")
io$trajectory <- "nmp"
io$trajectory_file <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacell_trajectory.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/pdf"); dir.create(io$outdir, showWarnings = F)

###################
## Load metadata ##
###################

metacell_metadata.dt <- fread(io$metacell_metadata)

##############
## Load RNA ##
##############

sce <- readRDS(io$metacell_sce)

#####################
## Load trajectory ##
#####################

# trajectory.dt <- fread(io$atlas_trajectory) %>% setnames(c("cell","V1","V2"))
trajectory.dt <- fread(io$trajectory_file) %>% 
  setnames(c("metacell","V1","V2")) %>% 
  merge(metacell_metadata.dt,by="metacell")

metacells <- intersect(colnames(sce), trajectory.dt$metacell)
trajectory.dt <- trajectory.dt[metacell%in%metacells]
sce <- sce[,colnames(sce)%in%metacells]

###############
## Load ATAC ##
###############

io$archR.peakMatrix.metacells <- file.path(io$basedir,"results/atac/archR/metacells/trajectories/nmp/PeakMatrix/PeakMatrix_summarized_experiment_metacells.rds")
atac_peakMatrix_metacells.se <- readRDS(io$archR.peakMatrix.metacells)
metacells.atac <- intersect(trajectory.dt$metacell,colnames(atac_peakMatrix_metacells.se))
atac_peakMatrix_metacells.se <- atac_peakMatrix_metacells.se[,metacells.atac]
assayNames(atac_peakMatrix_metacells.se) <- "counts"
assay(atac_peakMatrix_metacells.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_metacells.se),2,colSums(assay(atac_peakMatrix_metacells.se),na.rm=T),"/"))+1)

###########################################
## Plot trajectory coloured by cell type ##
###########################################

to.plot <- trajectory.dt

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  geom_point(aes(fill=celltype), size=2.5, shape=21, stroke=0.25) +
  # viridis::scale_fill_viridis() +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none"
  )

pdf(file.path(io$outdir,"metacell_trajectory_nmp_test.pdf"), width=6, height=4)
print(p)
dev.off()

#########################################################
## Plot trajectory coloured by one cell type at a time ##
#########################################################

to.plot <- trajectory.dt

celltypes.to.plot <- unique(to.plot$celltype)
for (i in celltypes.to.plot) {

  opts$colors <- opts$celltype.colors[celltypes.to.plot]
  opts$colors[names(opts$colors)!=i] <- "gray70"
  
  opts$alpha = rep(1,length(celltypes.to.plot))
  names(opts$alpha) <- celltypes.to.plot
  opts$alpha[names(opts$alpha)!=i] <- 0.4
  
  opts$size = rep(3,length(celltypes.to.plot))
  names(opts$size) <- celltypes.to.plot
  opts$size[names(opts$size)!=i] <- 1.25
  
  p <- ggplot(to.plot, aes(x=V1, y=V2)) +
    geom_point(aes(fill=celltype, size=celltype, alpha=celltype), size=2.5, shape=21, stroke=0.25) +
    # viridis::scale_fill_viridis() +
    scale_fill_manual(values=opts$colors) +
    scale_size_manual(values=opts$size) +
    scale_alpha_manual(values=opts$alpha) +
    labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position="none"
    )
  
  pdf(file.path(io$outdir,sprintf("metacell_trajectory_nmp_coloured_by_%s.pdf",i)), width=3.5, height=3.5)
  print(p)
  dev.off()
}


##############################
## Plot expression of genes ##
##############################
 
genes.to.plot <- c("Fgf4","Tbx6","Mesp1")

for (i in genes.to.plot) {
  rna.dt <- as.matrix(logcounts(sce[i,])) %>% 
    as.data.table(keep.rownames = "gene") %>% 
    melt(id.vars = "gene", variable.name = "metacell", value.name = "expr")
  
  to.plot <- trajectory.dt %>% merge(rna.dt, by="metacell", allow.cartesian=TRUE)
  to.plot[expr<=2.5,expr:=2.5]
  hist(to.plot$expr)
  
  p <- ggplot(to.plot, aes(x=V1, y=V2)) +
    geom_point(aes(fill=expr), size=2.5, shape=21, stroke=0.25) +
    # facet_wrap(~gene) +
    scale_fill_gradient(low = "gray95", high = "darkgreen") +
    labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position = "none"
    )
  
  pdf(file.path(io$outdir,sprintf("metacell_trajectory_nmp_coloured_by_%s.pdf",i)), width=5, height=3.5)
  print(p)
  dev.off()
}


##############################
## Plot interaction of genes #
##############################

geneA <- "T"
geneB <- "Sox2"

tmp <- data.table(
  metacell = colnames(sce),
  expr = logcounts(sce[geneA,])[1,]*logcounts(sce[geneB,])[1,] %>% minmax.normalisation
)

to.plot <- trajectory.dt %>% merge(tmp, by="metacell")

to.plot[expr<=1,expr:=0]
to.plot[expr>=6,expr:=6]

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  geom_point(aes(fill=expr), size=2.5, shape=21, stroke=0.25) +
  scale_fill_gradient(low = "gray80", high = "darkgreen") +
  labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none"
  )

pdf(file.path(io$outdir,"metacell_trajectory_nmp_T_Sox2_interaction.pdf"), width=5, height=3.5)
print(p)
dev.off()

#################################
## Plot accessibility of peaks ##
#################################

peaks.to.plot <- c("chr7:126785067-126785667", "chr7:79789147-79789747", "chr7:144884955-144885555")

for (i in peaks.to.plot) {
  
  acc.dt <- data.table(
    acc = assay(atac_peakMatrix_metacells.se,"logcounts")[i,],
    metacell = colnames(atac_peakMatrix_metacells.se),
    peak = i
  )
  
  to.plot <- trajectory.dt %>% merge(acc.dt, by="metacell")
  
  to.plot[acc<=1,acc:=1] # cap values for better viz
  # hist(to.plot$acc)
  
  p <- ggplot(to.plot, aes(x=V1, y=V2)) +
    geom_point(aes(fill=acc), size=2.5, shape=21, stroke=0.25) +
    # facet_wrap(~gene) +
    scale_fill_gradient(low = "gray95", high = "blue") +
    labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
    theme_classic() +
    ggplot_theme_NoAxes() +
    theme(
      legend.position = "none"
    )
  
  pdf(file.path(io$outdir,sprintf("metacell_trajectory_nmp_coloured_by_acc_%s.pdf",gsub(":","-",i))), width=5, height=3.5)
  print(p)
  dev.off()
}

#########################################
## Boxplots of chromatin accessibility ##
#########################################

peaks.to.plot <- c("chr7:144884955-144885555","chr7:79789147-79789747","chr7:126785067-126785667")

# i <- "chr7:144884955-144885555"
for (i in peaks.to.plot) {
  
  to.plot <- data.table(
    acc = assay(atac_peakMatrix_metacells.se,"logcounts")[i,],
    metacell = colnames(atac_peakMatrix_metacells.se)
  ) %>% merge(trajectory.dt,by="metacell")
  
  to.plot[celltype=="Caudal_Mesoderm",celltype:="Somitic_mesoderm"]
  celltype.order <- c("Spinal_cord","NMP","Somitic_mesoderm")
  to.plot[,celltype:=factor(celltype, levels=celltype.order)]
  
  my_comparisons <- list( c("Spinal_cord", "NMP"), c("NMP", "Somitic_mesoderm"))
  
  p <- ggboxplot(to.plot, x="celltype", y="acc", fill="celltype", outlier.shape=NA) +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), comparisons = my_comparisons, method="t.test") +
    stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width=0.75)) +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="", y="Chromatin accessibility (log normalised counts)") +
    # geom_violin(aes(fill=celltype)) +
    theme_classic() +
    theme(
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black"),
      legend.title = element_blank(),
      legend.position = "none"
    )
  
  pdf(file.path(io$outdir,sprintf("metacell_trajectory_nmp_boxplots_acc_%s.pdf",gsub(":","-",i))), width=5.5, height=6)
  print(p)
  dev.off()
}


