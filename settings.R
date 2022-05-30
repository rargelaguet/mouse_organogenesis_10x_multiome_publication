suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation_multiome_10x"
  io$atlas.basedir <- "/Users/ricard/data/gastrulation10x"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (Sys.info()[['nodename']]=="BI2404M") {
  io$basedir <- "/Users/argelagr/data/gastrulation_multiome_10x/test"
  io$atlas.basedir <- "/Users/argelagr/data/gastrulation10x"
  io$gene_metadata <- "/Users/argelagr/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$TFs_file <- "/Users/argelagr/data/mm10_regulation/TFs/TFs.txt"
} else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
  if (grepl("Clark", Sys.info()['effective_user'])) {
    io$basedir       <- "/bi/scratch/Stephen_Clark/multiome/resilio"
    io$rawdata       <- "/bi/scratch/Stephen_Clark/multiome/raw"
    io$gene_metadata <- "/bi/scratch/Stephen_Clark/annotations/Mmusculus_genes_BioMart.87.txt"
  } else if (grepl("argelag", Sys.info()['effective_user'])) {
    io$basedir <-'/bi/group/reik/ricard/data/gastrulation_multiome_10x/test'
    io$archR.directory <- file.path(io$basedir,"processed/atac/archR")
    io$gene_metadata <- "/bi/group/reik/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
    io$atlas.basedir <- "/bi/group/reik/ricard/data/pijuansala2019_gastrulation10x"
    io$TFs_file <- "/bi/group/reik/ricard/data/mm10_regulation/TFs/TFs.txt" 
  }
} else if (grepl("bi2228m", Sys.info()['nodename'])) {
  io$basedir      <- "/Users/clarks/data/multiome/"
  io$rawdata      <- "/Users/clarks/data/multiome_raw_sub/"
} else if(grepl('LAPTOP',Sys.info()['nodename'])){
  io$basedir      <- 'D:/OneDrive - zju.edu.cn/research/bioin/lab/Ricard/data/gastrulation_multiome_10x'
} else if(grepl('Workstation',Sys.info()['nodename'])){
  io$basedir      <-'/home/lijingyu/gastrulation/data/gastrulation_multiome_10x'
  io$gene_metadata <- "/home/lijingyu/gastrulation/data/gastrulation_multiome_10x/Mmusculus_genes_BioMart.87.txt"
  io$atlas.basedir <- "/home/lijingyu/gastrulation/data/gastrulation10x"
} else {
  stop("Computer not recognised")
}

io$metadata <- file.path(io$basedir,"sample_metadata.txt.gz")

# TFs
# io$TFs <- file.path(io$basedir,"results/TFs.txt")

# RNA
io$rna.anndata <- file.path(io$basedir,"processed/rna/anndata.h5ad")
io$rna.seurat <- file.path(io$basedir,"processed/rna/seurat.rds")
io$rna.sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
io$rna.tfs.sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment_TFs.rds")
io$rna.pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
io$rna.tfs.pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_TFs_pseudobulk.rds")
io$rna.metacells.sce <- file.path(io$basedir,"results/rna/metacells/all_cells/SingleCellExperiment_metacells.rds")

# Differential RNA expression and marker genes
# io$rna.differential <- file.path(io$basedir,"results/rna/differential")
io$rna.celltype.marker_genes.all <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_genes_all.txt.gz")
io$rna.celltype.marker_genes.filt <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_genes_filtered.txt.gz")
io$rna.celltype.marker_tfs.all <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_tfs_all.txt.gz")
io$rna.celltype.marker_tfs.filt <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_tfs_filtered.txt.gz")

# RNA atlas (Pijuan-Sala2019)
io$rna.atlas.metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
io$rna.atlas.sce <- file.path(io$atlas.basedir,"processed/SingleCellExperiment.rds")
io$rna.atlas.marker_genes.up <- file.path(io$atlas.basedir,"results/differential/celltypes/genes/all_stages/marker_genes/marker_genes_upregulated_filtered.txt.gz")
io$rna.atlas.marker_genes.all <- file.path(io$atlas.basedir,"results/differential/celltypes/genes/all_stages/marker_genes/marker_genes_upregulated_all.txt.gz")
io$rna.atlas.marker_TFs.up <- file.path(io$atlas.basedir,"results/differential/celltypes/TFs/TF_markers/marker_TFs_upregulated_filtered.txt.gz")
io$rna.atlas.marker_TFs.all <- file.path(io$atlas.basedir,"results/differential/celltypes/TFs/TF_markers/marker_TFs_upregulated_all.txt.gz")
io$rna.atlas.differential <- file.path(io$atlas.basedir,"results/differential")
# io$rna.atlas.average_expression_per_celltype <- file.path(io$atlas.basedir,"results/marker_genes/all_stages/avg_expr_per_celltype_and_gene.txt.gz")
io$rna.atlas.sce.pseudobulk <- file.path(io$atlas.basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk.rds")
io$rna.atlas.celltype_proportions <- file.path(io$atlas.basedir,"results/celltype_proportions/celltype_proportions.txt.gz")

# ATAC
# io$signac.directory <- file.path(io$basedir,"processed/atac/signac")
# io$signac <- file.path(io$signac.directory,"/signac.rds")
# io$cistopic.model <- file.path(io$basedir,"results/atac/cistopic/cistopic_warpLDA_bestModel.rds")
io$archR.directory <- file.path(io$basedir,"processed/atac/archR")
io$motifmatcher_scores.se <- file.path(io$archR.directory,"Annotations/CISBP-Scores.rds")
io$motifmatcher_positions.se <-file.path(io$archR.directory,"Annotations/CISBP-Positions.rds")
io$archR.projectMetadata <- file.path(io$archR.directory,"projectMetadata.rds")
io$archR.peakSet.granges <- file.path(io$archR.directory,"PeakSet.rds")
io$archR.peakSet.stats <- file.path(io$basedir,"results/atac/archR/feature_stats/PeakMatrix/PeakMatrix_celltype_stats.txt.gz")
io$archR.bgdPeaks <- file.path(io$archR.directory,"Background-Peaks.rds")
io$archR.peakSet.bed <- file.path(io$archR.directory,"PeakCalls/peaks_archR_macs2.bed.gz")
io$archR.peakMatrix.cells <- file.path(io$archR.directory,"Matrices/PeakMatrix_summarized_experiment.rds")
io$archR.peakMatrix.metacells <- file.path(io$basedir,"results/atac/archR/metacells/all_cells/PeakMatrix/PeakMatrix_summarized_experiment_metacells.rds")
io$archR.peakMatrix.pseudobulk <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype/PeakMatrix/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$archR.peak.metadata <- file.path(io$archR.directory,"PeakCalls/peak_metadata.tsv.gz")
io$archR.peak.stats <- file.path(io$basedir,"results/atac/archR/feature_stats/PeakMatrix_celltype.mapped_stats.txt.gz")
io$archR.peak2gene.all <- file.path(io$basedir,"results/atac/archR/peak_calling/peaks2genes/peaks2genes_all.txt.gz")
io$archR.peak2gene.nearest <- file.path(io$basedir,"results/atac/archR/peak_calling/peaks2genes/peaks2genes_nearest.txt.gz")
io$archR.GeneScoreMatrix.cells <- file.path(io$basedir,"processed/atac/archR/Matrices/GeneScoreMatrix_TSS_summarized_experiment.rds")
io$archR.GeneScoreMatrix.pseudobulk <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype/GeneScoreMatrix_TSS/pseudobulk_GeneScoreMatrix_TSS_summarized_experiment.rds")
io$archR.chromvar.cells <- file.path(io$basedir,"results/atac/archR/chromvar/cells/chromVAR_deviations_CISBP_archr.rds")
io$archR.chromvar.pseudobulk <- file.path(io$basedir,"results/atac/archR/chromvar/pseudobulk/celltype/chromVAR_deviations_CISBP_pseudobulk_archr.rds")
io$archR.feature_stats <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_CISBP_archr.rds")
io$archR.motif2gene <- file.path(io$basedir,"processed/atac/archR/Annotations/CISBP_motif2gene.txt.gz")
# io$archR.chromvar_chip.cells <- file.path(io$basedir,"results/atac/archR/chromvar/cells/celltype/chromVAR_deviations_CISBP_pseudobulk_archr.rds")

# io$archR.pseudobulk.genes.se <- file.path(io$basedir,"processed/atac/archR/pseudobulk/celltype.predicted/pseudobulk_GeneScoreMatrix_summarized_experiment.rds")

# Differential ATAC and marker peaks
io$atac.markers_peaks.pseudobulk.all <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype/PeakMatrix/parsed/markers_all.txt.gz")
io$atac.markers_peaks.pseudobulk.filt <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype/PeakMatrix/parsed/markers_filt.txt.gz")
# io$archR.markers_chromvar <- file.path(io$basedir,"results/atac/archR/chromvar/differential/markers/marker_TFs.txt.gz")

# ATAC public data (Pijuan-Sala2020)
io$pijuansala.basedir <- file.path(io$basedir,"public_datasets/Pijuan-Sala_2020")
io$pijuansala.archR.directory <- file.path(io$pijuansala.basedir,"data/processed/archR")

# paga
io$paga.connectivity <- file.path(io$atlas.basedir,"results/paga/paga_connectivity.csv")
io$paga.coordinates <- file.path(io$atlas.basedir,"results/paga/paga_coordinates.csv")

# Virtual ChIP-seq
# io$tf2peak_cor.se <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
# io$tf.activities <- file.path(io$atlas.basedir,"results/rna_atac/rna_vs_chromvar/TF_activities/tf_activities.txt.gz")
io$virtual_chip.dir <- file.path(io$basedir,"results/rna_atac/virtual_chipseq")
io$virtual_chip.mtx <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/virtual_chip_mtx.rds")
io$archR.chromvar_chip.pseudobulk <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_CISBP_archr.rds")

# Coexpression
# io$tf2tf_cor.mtx <- file.path(io$basedir,"results/rna/coexpression/correlation_matrix_tf2tf.rds")
# io$tf2gene_cor.mtx <- file.path(io$basedir,"results/rna/coexpression/correlation_matrix_tf2gene.rds")

# Dimensionality reduction
io$rna.dimred.pca <- file.path(io$basedir,"results/rna/dimensionality_reduction/sce/batch_correction_sample_remove_ExE_cells_False/pca_features2500_pcs50.txt.gz")
io$atac.dimred.lsi <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/cells/PeakMatrix/remove_ExE_cells_False/batch_correction_None/lsi_nfeatures25000_ndims50.txt.gz")
# io$rna.dimred.umap <- file.path(io$basedir,"results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_features2500_pcs30_neigh25_dist0.3.txt.gz")
# io$atac.dimred.umap <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims50_neigh30_dist0.45.txt.gz")

#############
## Options ##
#############

opts <- list()

opts$stages <- c(
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.5",
  "E8.75"
)

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E7.75_rep1",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2",
  "E8.5_CRISPR_T_WT",
  "E8.5_CRISPR_T_KO",
  "E8.75_rep1",
  "E8.75_rep2"
)

# opts$rename.samples <- c(
#   "E8.75_rep1" = "e8_75_1_L002",
#   "E8.75_rep2" = "e8_75_2_L002",
#   "E8.5_rep1" = "multiome1",
#   "E8.5_rep2" = "multiome2",
#   "E8.0_rep1" = "E8_0_rep1_multiome",
#   "E8.0_rep2" = "E8_0_rep2_multiome",
#   "E7.5_rep1" = "rep1_L001_multiome",
#   "E7.5_rep2" = "rep2_L002_multiome"
#   # E7.75_rep1
#   # E8.5_CRISPR_T_WT
#   # E8.5_CRISPR_T_KO
# )

opts$sample2stage <- c(
  "E7.5_rep1" = "E7.5",
  "E7.5_rep2" = "E7.5",
  "E7.75_rep1" = "E7.75",
  "E8.0_rep1" = "E8.0",
  "E8.0_rep2" = "E8.0",
  "E8.5_rep1" = "E8.5",
  "E8.5_rep2" = "E8.5",
  "E8.75_rep1" = "E8.75",
  "E8.75_rep2" = "E8.75",
  "E8.5_CRISPR_T_KO" = "E8.5",
  "E8.5_CRISPR_T_WT" = "E8.5"
)


opts$celltypes <- c(
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

opts$stage.colors = c(
  "E7.5" = "#edf8b1",
  "E7.75" = "#c7e9b4",
  "E8.0" = "#7fcdbb",
  "E8.25" = "#41b6c4",
  "E8.5" = "#1d91c0",
  "E8.75" = "#225ea8"
)
# opts$stage.colors <- viridis::viridis(n=length(opts$stages))
# names(opts$stage.colors) <- rev(opts$stages)

opts$celltype.colors = c(
  "Epiblast" = "#635547",
  "Primitive_Streak" = "#DABE99",
  "Caudal_epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  "Anterior_Primitive_Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def._endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent_mesoderm" = "#C594BF",
  "Mixed_mesoderm" = "#DFCDE4",
  "Intermediate_mesoderm" = "#139992",
  "Caudal_Mesoderm" = "#3F84AA",
  "Paraxial_mesoderm" = "#8DB5CE",
  "Somitic_mesoderm" = "#005579",
  "Pharyngeal_mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE_mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial_progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  "Blood_progenitors" = "#c9a997",
  "Blood_progenitors_1" = "#f9decf",
  "Blood_progenitors_2" = "#c9a997",
  "Erythroid" = "#EF4E22",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid3" = "#EF4E22",
  "NMP" = "#8EC792",
  "Neurectoderm" = "#65A83E",
  "Rostral_neurectoderm" = "#65A83E",
  "Caudal_neurectoderm" = "#354E23",
  "Neural_crest" = "#C3C388",
  "Forebrain_Midbrain_Hindbrain" = "#647a4f",
  "Spinal_cord" = "#CDE088",
  "Surface_ectoderm" = "#f7f79e",
  "Visceral_endoderm" = "#F6BFCB",
  "ExE_endoderm" = "#7F6874",
  "ExE_ectoderm" = "#989898",
  "Parietal_endoderm" = "#1A1A1A"
)

opts$chr <- paste0("chr",c(1:19,"X","Y"))

###################
## Load metadata ##
###################

# sample_metadata <- fread(io$metadata)
# .[pass_QC==T] %>% 
# .[batch%in%opts$batches] %>%
# .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ","_")] %>%
# .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/","_")] %>%
# .[,celltype.mapped:=factor(celltype.mapped, levels=names(opts$celltype.colors))]


###################
## Edit metadata ##
###################

# io$metadata <- file.path(io$basedir,"results/atac/archR/qc/sample_metadata_after_qc.txt.gz")
# io$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# sample_metadata <- fread(io$metadata)
# sample_metadata[pass_rnaQC==FALSE & nFrags_atac<=1e4,pass_atacQC:=FALSE]
# sample_metadata[sample=="E7.75_rep1" & ribosomal_percent_RNA>=8,pass_rnaQC:=FALSE]
# sample_metadata[sample=="E7.75_rep1" & ribosomal_percent_RNA>=4.5 & nFeature_RNA<=3500,pass_rnaQC:=FALSE]
# fwrite(sample_metadata, io$metadata, sep="\t", quote=F, na="NA")


