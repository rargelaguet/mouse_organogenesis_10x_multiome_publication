# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O metadata
io$cell_metadata <- io$metadata
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacells_metadata.txt.gz")

# I/O TF expr vs peak acc
io$metacells_tf_expr_vs_peak_acc <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/trajectories/nmp/TFexpr_vs_peakAcc/CISBP_cor_TFexpr_vs_peakAcc.rds")

# I/O RNA
io$cells_sce <- io$rna.sce
io$metacells_sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")

# I/O ATAC
io$cells_atac_peak_matrix <- file.path(io$basedir,"processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
io$metacells_atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/metacells/trajectories/nmp/PeakMatrix_summarized_experiment_metacells.rds")

# I/O Output directory
io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/trajectories/nmp"); dir.create(io$outdir, showWarnings = F)

# Options
opts$samples <- c("E8.5_CRISPR_T_KO", "E8.5_CRISPR_T_WT")
opts$celltypes <- c("Caudal_Mesoderm", "Somitic_mesoderm", "NMP", "Spinal_cord")

####################
## Load metadata  ##
####################

metacell_metadata.dt <- fread(io$metacell_metadata) %>%
  .[sample%in%opts$samples & celltype%in%opts$celltypes]

cell_metadata.dt <- fread(io$cell_metadata) %>% 
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype%in%opts$celltypes]

###########################################################
## Load TF expression vs peak accessibility correlations ##
###########################################################

tf_expr_vs_peak_acc_metacells.se <- readRDS(io$metacells_tf_expr_vs_peak_acc)

opts$TFs <- colnames(tf_expr_vs_peak_acc_metacells.se)

########################
## Load metacell data ##
########################

# Load RNA SingleCellExperiment
rna_metacells.sce <- readRDS(io$metacells_sce)[,metacell_metadata.dt$metacell]
rna_metacells.sce <- rna_metacells.sce[toupper(rownames(rna_metacells.sce))%in%opts$TFs]
rownames(rna_metacells.sce) <- toupper(rownames(rna_metacells.sce))

# Load ATAC SummarizedExperiment
atac_peakMatrix_metacells.se <- readRDS(io$metacells_atac_peak_matrix)

# Normalise atac data
# assay(atac_peakMatrix_metacells.se)

# Make sure that samples are consistent
samples <- intersect(colnames(rna_metacells.sce),colnames(atac_peakMatrix_metacells.se))
rna_metacells.sce <- rna_metacells.sce[,samples]
atac_peakMatrix_metacells.se <- atac_peakMatrix_metacells.se[,samples]
metacell_metadata.dt <- metacell_metadata.dt[metacell%in%samples] %>% setkey(metacell) %>% .[samples]

####################
## Load cell data ##
####################

# Load RNA SingleCellExperiment
rna_cells.sce <- load_SingleCellExperiment(file=io$cells_sce, cells=cell_metadata.dt$cell, remove_non_expressed_genes = T, normalise = T)
colData(rna_cells.sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame
rna_cells.sce <- rna_cells.sce[toupper(rownames(rna_cells.sce))%in%opts$TFs]
rownames(rna_cells.sce) <- toupper(rownames(rna_cells.sce))

# Load ATAC SummarizedExperiment
atac_peakMatrix_cells.se <- readRDS(io$cells_atac_peak_matrix)[,cell_metadata.dt$cell]

##########
## Plot ##
##########

stopifnot(colnames(rna_metacells.sce)==colnames(atac_peakMatrix_metacells.se))

genes.to.plot <- c("T")

# i <- "T"
for (i in genes.to.plot) {
  
  tmp <- data.table(
    gene = i,
    peak = rownames(tf_expr_vs_peak_acc_metacells.se),
    metacells_cor = assay(tf_expr_vs_peak_acc_metacells.se[,i],"cor")[,1]
  )
  
  peaks.to.plot <- sort(assay(tf_expr_vs_peak_acc_metacells.se[,i],"cor")[,1]) %>% tail(n=5) %>% names
  
  # j <- "chr7:79789167-79789767"# peaks.to.plot[1]
  for (j in peaks.to.plot) {

    cells_expr_acc.dt <- data.table(
      gene = i,
      cell = colnames(rna_cells.sce),
      celltype = cell_metadata.dt$celltype,
      expr = logcounts(rna_cells.sce[i,])[1,],
      acc = assay(atac_peakMatrix_cells.se)[j,]
    )

    p.cells <- ggplot(cells_expr_acc.dt, aes(x=acc, y=expr, color=celltype)) +
      geom_point(size=0.5) +
      stat_smooth(method="lm") +
      scale_fill_manual(values=opts$celltype.colors) +
      labs(x="Chromatin accessibility", y="Gene expression") +
      theme_classic() +
      theme(
        axis.text = element_text(colour="black",size=rel(0.8)),
        legend.position = "none"
      )
    
    metacells_expr_acc.dt <- data.table(
      gene = i,
      cell = colnames(rna_metacells.sce),
      celltype = metacell_metadata.dt$celltype,
      expr = logcounts(rna_metacells.sce[i,])[1,],
      acc = assay(atac_peakMatrix_metacells.se)[j,]
    )
    
    p.metacells <- ggplot(metacells_expr_acc.dt, aes(x=acc, y=expr)) +
      geom_point(aes(fill=celltype), shape=21) +
      scale_fill_manual(values=opts$celltype.colors) +
      labs(x="Chromatin accessibility (metacells)", y="Gene expression (metacells)") +
      stat_smooth(method="lm", color="black") +
      stat_cor(method = "pearson") +
      theme_classic() +
      theme(
        axis.text = element_text(colour="black",size=rel(0.8)),
        legend.position = "none"
      )
    
    p <- cowplot::plot_grid(plotlist=list(p.cells,p.metacells), ncol=2, rel_widths = c(0.5,0.5))
    
    pdf(file.path(io$outdir,sprintf("%s_cell_metacell_tf_expr_vs_peak_acc.pdf",i)), width=10, height=10)
    print(p)
    dev.off()
    
  }
}


##########
## TEST ##
##########

j <- "chr7:79789167-79789767"
max(assay(atac_peakMatrix_metacells.se))
