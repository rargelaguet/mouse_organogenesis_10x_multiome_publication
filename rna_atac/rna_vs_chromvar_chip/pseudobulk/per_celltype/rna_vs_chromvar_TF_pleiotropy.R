source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"
# opts$celltypes <- c(
#   "Epiblast",
#   "Primitive_Streak",
#   "Caudal_epiblast",
#   "PGC",
#   "Anterior_Primitive_Streak",
#   "Notochord",
#   "Def._endoderm",
#   "Gut",
#   "Nascent_mesoderm",
#   "Mixed_mesoderm",
#   "Intermediate_mesoderm",
#   "Caudal_Mesoderm",
#   "Paraxial_mesoderm",
#   "Somitic_mesoderm",
#   "Pharyngeal_mesoderm",
#   "Cardiomyocytes",
#   "Allantois",
#   "ExE_mesoderm",
#   "Mesenchyme",
#   "Haematoendothelial_progenitors",
#   "Endothelium",
#   "Blood_progenitors_1",
#   "Blood_progenitors_2",
#   "Erythroid1",
#   "Erythroid2",
#   "Erythroid3",
#   "NMP",
#   "Rostral_neurectoderm",
#   "Caudal_neurectoderm",
#   "Neural_crest",
#   "Forebrain_Midbrain_Hindbrain",
#   "Spinal_cord",
#   "Surface_ectoderm",
#   "Visceral_endoderm",
#   "ExE_endoderm",
#   "ExE_ectoderm",
#   "Parietal_endoderm"
# )

# I/O
io$rna_sce_pseudobulk <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_TFs_pseudobulk.rds")
io$atac_chromvar_chip_pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_archr.rds",opts$motif_annotation))
# io$atac_chromvar_chip_pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_chromvar.rds",opts$motif_annotation))
io$tf_markers_rna <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_tfs_all.txt.gz")
io$tf_markers_atac <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltypes/CISBP/parsed/markers_all.txt.gz")
io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_celltype/pleiotropy/%s",opts$motif_annotation)); dir.create(io$outdir, showWarnings=F, recursive = T)

###########################
## Load TF marker scores ##
###########################

tf_markers_rna.dt <- fread(io$tf_markers_rna) %>% 
  .[celltype%in%opts$celltypes] %>%
  setnames("score","tf_marker_score_rna")
tf_markers_atac.dt <- fread(io$tf_markers_atac) %>% 
  .[celltype%in%opts$celltypes] %>%
  setnames("score","tf_marker_score_atac")

tf_markers.dt <- merge(tf_markers_rna.dt, tf_markers_atac.dt, by=c("celltype","gene"))

tf_markers.dt[,tf_marker_score:=minmax.normalisation(tf_marker_score_rna*tf_marker_score_atac)] %>% setorder(-tf_marker_score)


#############################################
## Plot number of marker TFs per cell type ##
#############################################

celltypes.to.plot <- c("Somitic_mesoderm","Endothelium","Cardiomyocytes","Blood_progenitors_2","NMP","Endothelium","Spinal_cord")

seq.ranges <- seq(0.75,1,by=0.05)

to.plot <- seq.ranges %>% map(function(i) {
  tf_markers_rna.dt[tf_marker_score_rna>=i,.N,by="celltype"] %>% .[,min_score:=i]
}) %>% rbindlist %>% .[celltype%in%celltypes.to.plot]

p <- ggplot(to.plot, aes_string(x="min_score", y="N", color="celltype")) +
  geom_line(size=1.25) +
  labs(y="Number of TF markers", x="Minimum TF marker score (RNA)") +
  scale_color_manual(values=opts$celltype.colors[celltypes.to.plot]) +
  # coord_cartesian(xlim=c(0,0.9)) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black"),
    # legend.position = "right",
    legend.position = "none",
  )

pdf(file.path(io$outdir,"lineplot_number_TF_markers_per_celltype.pdf"), width=9.5, height=3)
print(p)
dev.off()


######################################
## Plot number of cell types per TF ##
######################################

to.plot <- tf_markers_rna.dt[tf_marker_score_rna>=0.75,.N,by="gene"] 

p <- ggplot(to.plot, aes_string(x="min_score", y="N", color="celltype")) +
  geom_line(size=1.25) +
  labs(y="Number of TF markers", x="Minimum TF marker score (RNA)") +
  scale_color_manual(values=opts$celltype.colors[celltypes.to.plot]) +
  # coord_cartesian(xlim=c(0,0.9)) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black"),
    # legend.position = "right",
    legend.position = "none",
  )

pdf(file.path(io$outdir,"lineplot_number_TF_markers_per_celltype.pdf"), width=9.5, height=3)
print(p)
dev.off()
