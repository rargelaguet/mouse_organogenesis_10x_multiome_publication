source(here::here("settings.R"))
source(here::here("utils.R"))

library(pheatmap)

#####################
## Define settings ##
#####################

opts$motif_annotation <- "CISBP"
io$rna_sce_pseudobulk <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_TFs_pseudobulk.rds")
io$atac_chromvar_chip_pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_archr.rds",opts$motif_annotation))
# io$atac_chromvar_chip_pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_chromvar.rds",opts$motif_annotation))
io$tf_markers_rna <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_tfs_all.txt.gz")
io$tf_markers_atac <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltypes/CISBP/parsed/markers_all.txt.gz")
io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_celltype/%s",opts$motif_annotation))

# I/O
dir.create(io$outdir, showWarnings=F, recursive = T)

# Options
# opts$celltypes <- c(
#   # "Epiblast",
#   # "Primitive_Streak",
#   # "Caudal_epiblast",
#   # "PGC",
#   # "Anterior_Primitive_Streak",
#   "Notochord",
#   "Def._endoderm",
#   "Gut",
#   "Nascent_mesoderm",
#   # "Mixed_mesoderm",
#   "Intermediate_mesoderm",
#   # "Caudal_Mesoderm",
#   "Paraxial_mesoderm",
#   "Somitic_mesoderm",
#   "Pharyngeal_mesoderm",
#   "Cardiomyocytes",
#   "Allantois",
#   # "ExE_mesoderm",
#   "Mesenchyme",
#   "Haematoendothelial_progenitors",
#   "Endothelium",
#   "Blood_progenitors_1",
#   "Blood_progenitors_2",
#   # "Erythroid1",
#   # "Erythroid2",
#   "Erythroid3",
#   "NMP",
#   "Rostral_neurectoderm",
#   # "Caudal_neurectoderm",
#   "Neural_crest",
#   "Forebrain_Midbrain_Hindbrain",
#   "Spinal_cord",
#   "Surface_ectoderm"
#   # "Visceral_endoderm",
#   # "ExE_endoderm",
#   # "ExE_ectoderm",
#   # "Parietal_endoderm"
# )

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load pseudobulk TF RNA expression
rna_pseudobulk.sce <- readRDS(io$rna_sce_pseudobulk)[,opts$celltypes]

# Load pseudobulk ATAC chromVAR-ChIP matrix
atac_pseudobulk_chromvar.se <- readRDS(io$atac_chromvar_chip_pseudobulk)[,opts$celltypes]

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

#############
## Heatmap ##
#############

# TFs.to.highlight <- tf_markers.dt %>% split(.$celltype) %>% map(function(x) head(x, n=2) %>% .$gene) %>% unlist %>% unique
# TFs.to.highlight <- c(TFs.to.highlight,c("TAL1","FOXA2","RFX4")) %>% unique
TFs.to.highlight <- c("TAL1","FOXA2","RFX4")

TFs.to.plot <- unique(tf_markers.dt[tf_marker_score>=0.75,gene])

to.plot <- tf_markers.dt %>% 
  .[gene%in%TFs.to.plot] %>%
  # .[tf_marker_score>=0.25] %>%
  dcast(celltype~gene, value.var="tf_marker_score", fill=0) %>%
  matrix.please

to.plot <- to.plot[opts$celltypes[opts$celltypes%in%rownames(to.plot)],]

rownames(to.plot) <- str_replace_all(rownames(to.plot),"_"," ")
colnames(to.plot)[!colnames(to.plot)%in%TFs.to.highlight] <- ""
# colnames(to.plot)[!seq(1,ncol(to.plot))%in%seq(1,ncol(to.plot),5)] <- ""

pheatmap(
  mat = to.plot, 
  cluster_cols = T, cluster_rows = F, 
  color = colorRampPalette(c("#F2F2F2", "#CD0000"))(100),
  fontsize_row = 7, fontsize_col = 5,
  treeheight_row = 0, treeheight_col = 0,
  legend = FALSE,
  # angle_col = 0,
  # scale = "column",
  filename = file.path(io$outdir,"heatmap_marker_TFs.pdf"),
  width = 7, height = 4
)

