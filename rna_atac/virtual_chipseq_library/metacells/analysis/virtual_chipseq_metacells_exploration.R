
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"
opts$trajectory <- "nmp"

# I/O
io$rna_metacells.sce <- file.path(io$basedir, 'results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds')
io$metacell_metadata <- file.path(io$basedir, 'results/atac/archR/metacells/trajectories/nmp/PeakMatrix/metacells_metadata.txt.gz')
io$archR.peakMatrix.metacells <- file.path(io$basedir,"results/atac/archR/metacells/trajectories/nmp/PeakMatrix/PeakMatrix_summarized_experiment_metacells.rds")
io$virtual_chip.mtx <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/%s/virtual_chip_matrix.rds",opts$motif_annotation))
io$trajectory <- "nmp"
io$trajectory_file <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacell_trajectory.txt.gz")
io$outdir <-  file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/%s/pdf",opts$trajectory)); dir.create(io$outdir, showWarnings = F)


if (io$trajectory=="nmp") {
  celltypes.to.plot <- c("Caudal_Mesoderm", "Somitic_mesoderm", "NMP", "Spinal_cord")
}

###################
## Load metadata ##
###################

metadata.dt <- fread(io$metacell_metadata) %>% .[celltype%in%celltypes.to.plot]

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(io$trajectory_file) %>% setnames(c("metacell","V1","V2"))

##################################
## Load virtual ChIP-seq matrix ##
##################################

virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)

##################################
## Load chromatin accessibility ##
##################################

atac_peakMatrix_metacells.se <- readRDS(io$archR.peakMatrix.metacells)

metacells <- intersect(trajectory.dt$metacell,colnames(atac_peakMatrix_metacells.se))

# Normalise ATAC data
assayNames(atac_peakMatrix_metacells.se) <- "counts"
assay(atac_peakMatrix_metacells.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix_metacells.se),2,colSums(assay(atac_peakMatrix_metacells.se),na.rm=T),"/"))+1)

##########
## Plot ##
##########

brachyury_binding_sites <- virtual_chip.mtx[,"T"][virtual_chip.mtx[,"T"]>=0.40] %>% sort

tmp <- assay(atac_peakMatrix_metacells.se,"logcounts")[names(brachyury_binding_sites),]

to.plot <- data.table(
  acc = colMeans(tmp),
  metacell = colnames(tmp)
) %>% merge(trajectory.dt,by="metacell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  geom_point(aes(fill=acc), size=2.5, shape=21, stroke=0.25) +
  # facet_wrap(~gene) +
  scale_fill_gradient(low = "gray95", high = "darkgreen") +
  labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.position = "right"
  )

pdf(file.path(io$outdir,sprintf("network_coloured_by_%s_expr.pdf",i)), width = 5, height = 5.5)
print(p)
dev.off()


