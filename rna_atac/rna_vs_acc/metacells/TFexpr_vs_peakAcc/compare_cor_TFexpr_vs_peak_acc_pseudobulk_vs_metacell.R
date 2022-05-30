# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O metadata
io$cell_metadata <- io$metadata
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/metacells_metadata.txt.gz")

# I/O TF expr vs peak acc
io$metacell_tf_expr_vs_peak_acc <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/CISBP_cor_TFexpr_vs_peakAcc.rds")
# io$metacells_tf_expr_vs_peak_acc <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/CISBP_cor_TFexpr_vs_peakAcc.rds")
io$metacells_tf_expr_vs_peak_acc <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/test/CISBP_cor_TFexpr_vs_peakAcc.rds")
io$pseudobulk_tf_expr_vs_peak_acc <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/CISBP_cor_TFexpr_vs_peakAcc.rds")

# I/O RNA
io$cells_sce <- io$rna.sce
io$metacells_sce <- file.path(io$basedir,"results/rna/metacells/SingleCellExperiment_metacells.rds")
io$pseudobulk_sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype.mapped/SingleCellExperiment_pseudobulk.rds")

# I/O ATAC
io$cells_atac_peak_matrix <- file.path(io$basedir,"processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
io$metacells_atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/metacells/PeakMatrix_summarized_experiment_metacells.rds")
io$pseudobulk_atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype.mapped/pseudobulk_PeakMatrix_summarized_experiment.rds")

# I/O Output directory
# io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc"); dir.create(io$outdir, showWarnings = F)
io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/pdf")

dir.create(io$outdir, showWarnings = F)
dir.create(file.path(io$outdir,"individual_examples"), showWarnings = F)

# Options
# opts$samples <- c("E8.5_CRISPR_T_KO", "E8.5_CRISPR_T_WT")
# opts$celltypes <- c("Caudal_Mesoderm", "Somitic_mesoderm", "NMP", "Spinal_cord")

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

tf_expr_vs_peak_acc_metacells.se <- readRDS(io$metacell_tf_expr_vs_peak_acc)
tf_expr_vs_peak_acc_pseudobulk.se <- readRDS(io$pseudobulk_tf_expr_vs_peak_acc)

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

# Normalise ATAC data
assay(atac_peakMatrix_metacells.se) <- log2(1e6*(sweep(assay(atac_peakMatrix_metacells.se),2,colSums(assay(atac_peakMatrix_metacells.se)),"/"))+0.5)
# assay(atac_peakMatrix_metacells.se)[1:5,1:5]

# Make sure that samples are consistent
samples <- intersect(colnames(rna_metacells.sce),colnames(atac_peakMatrix_metacells.se))
rna_metacells.sce <- rna_metacells.sce[,samples]
atac_peakMatrix_metacells.se <- atac_peakMatrix_metacells.se[,samples]
metacell_metadata.dt <- metacell_metadata.dt[metacell%in%samples] %>% setkey(metacell) %>% .[samples]

##########################
## Load pseudobulk data ##
##########################

# Load RNA SingleCellExperiment
rna_pseudobulk.sce <- readRDS(io$pseudobulk_sce)
rna_pseudobulk.sce <- rna_pseudobulk.sce[toupper(rownames(rna_pseudobulk.sce))%in%opts$TFs]
rownames(rna_pseudobulk.sce) <- toupper(rownames(rna_pseudobulk.sce))

# Load ATAC SummarizedExperiment
atac_peakMatrix_pseudobulk.se <- readRDS(io$pseudobulk_atac_peak_matrix)

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

###########
## Parse ##
###########

# # Subset celltypes
# cell.celltypes <- cell_metadata.dt[,.N,by=c("celltype")] %>% .[N>=5,celltype]
# metacell.celltypes <- metacell_metadata.dt[,.N,by=c("celltype")] %>% .[N>=5,celltype]
# celltypes.to.plot <- opts$celltypes[opts$celltypes%in%intersect(cell.celltypes,metacell.celltypes)]
# 
# # Subset SingleCellExperiment
# cells.sce <- cells.sce[,cells.sce$celltype%in%celltypes.to.plot]
# metacells.sce <- metacells.sce[,metacells.sce$celltype%in%celltypes.to.plot]
# pseudobulk.sce <- pseudobulk.sce[,colnames(pseudobulk.sce)%in%celltypes.to.plot]

################################################
## Plot histogram of correlation coefficients ##
################################################

tf_expr_vs_peak_acc_pseudobulk.cor <- as.numeric(assay(tf_expr_vs_peak_acc_pseudobulk.se,"cor"))
tf_expr_vs_peak_acc_pseudobulk.cor <- tf_expr_vs_peak_acc_pseudobulk.cor[tf_expr_vs_peak_acc_pseudobulk.cor!=0]

tf_expr_vs_peak_acc_metacells.cor <- as.numeric(assay(tf_expr_vs_peak_acc_metacells.se,"cor"))
tf_expr_vs_peak_acc_metacells.cor <- tf_expr_vs_peak_acc_metacells.cor[tf_expr_vs_peak_acc_metacells.cor!=0]

tmp <- sample.int(length(tf_expr_vs_peak_acc_metacells.cor), size = 1e4)
to.plot <- rbind(
  data.table(class = "metacell", cor = tf_expr_vs_peak_acc_metacells.cor[tmp]),
  data.table(class = "pseudobulk", cor = tf_expr_vs_peak_acc_pseudobulk.cor[tmp])
)

p <- gghistogram(to.plot, x="cor", fill="class", bins=100) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(file.path(io$outdir,"tf_expr_vs_peak_acc_cor_histogram_metacell_vs_pseudobulk.pdf"), width=5, height=4)
print(p)
dev.off()

##########################################
## Plot individual TF2peak correlations ##
##########################################

stopifnot(colnames(rna_metacells.sce)==colnames(atac_peakMatrix_metacells.se))

genes.to.plot <- c("GATA1","T","FOXA2")

# i <- "GATA1"
for (i in genes.to.plot) {
  
  tmp <- data.table(
    gene = i,
    peak = rownames(tf_expr_vs_peak_acc_pseudobulk.se),
    pseudobulk_cor = assay(tf_expr_vs_peak_acc_pseudobulk.se[,i],"cor")[,1],
    metacells_cor = assay(tf_expr_vs_peak_acc_metacells.se[,i],"cor")[,1]
  )
  
  peaks.to.plot <- sort(assay(tf_expr_vs_peak_acc_pseudobulk.se[,i],"cor")[,1]) %>% tail(n=5) %>% names
  
  # j <- "chr1:71769451-71770051"# peaks.to.plot[1]
  for (j in peaks.to.plot) {

    # verify correlation
    # assay(tf_expr_vs_peak_acc_pseudobulk.se[j,i])
    # assay(tf_expr_vs_peak_acc_metacells.se[j,i])
    # cor(logcounts(rna_metacells.sce[i,])[1,], assay(atac_peakMatrix_metacells.se)[j,])
    # cor(logcounts(rna_pseudobulk.sce[i,])[1,], assay(atac_peakMatrix_pseudobulk.se)[j,])
    
    # cells_expr_acc.dt <- data.table(
    #   gene = i,
    #   cell = colnames(rna_cells.sce),
    #   celltype = cell_metadata.dt$celltype,
    #   expr = logcounts(rna_cells.sce[i,])[1,],
    #   acc = assay(atac_peakMatrix_cells.se)[j,]
    # )
    # 
    # p.cells <- ggplot(cells_expr_acc.dt, aes(x=acc, y=expr, color=celltype)) +
    #   geom_point(size=0.5) +
    #   stat_smooth(method="lm") +
    #   scale_fill_manual(values=opts$celltype.colors) +
    #   labs(x="Chromatin accessibility", y="Gene expression") +
    #   theme_classic() +
    #   theme(
    #     axis.text = element_text(colour="black",size=rel(0.8)),
    #     legend.position = "none"
    #   )
    
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
    
    pseudobulk_expr_acc.dt <- data.table(
      gene = i,
      celltype = colnames(rna_pseudobulk.sce),
      expr = logcounts(rna_pseudobulk.sce[i,])[1,],
      acc = assay(atac_peakMatrix_pseudobulk.se)[j,]
    )
    
    p.pseudobulk <- ggplot(pseudobulk_expr_acc.dt, aes(x=acc, y=expr)) +
      geom_point(aes(fill=celltype), shape=21, size=4) +
      scale_fill_manual(values=opts$celltype.colors) +
      stat_smooth(method="lm", color="black") +
      stat_cor(method = "pearson") +
      labs(x="Chromatin accessibility (pseudobulk)", y="Gene expression (pseudobulk)") +
      theme_classic() +
      theme(
        axis.text = element_text(colour="black",size=rel(0.8)),
        legend.position = "none"
      )
    
    p <- cowplot::plot_grid(plotlist=list(p.metacells,p.pseudobulk), ncol=2, rel_widths = c(0.5,0.5))
    
    pdf(file.path(io$outdir,sprintf("individual_examples/%s_%s_cell_metacell_pseudobulk_expr.pdf",i,gsub("[:-]","_",j))), width=10, height=5)
    print(p)
    dev.off()
    
  }
}


##########
## TEST ##
##########

# j <- "chr7:79789167-79789767"
# max(assay(atac_peakMatrix_metacells.se))
# max(assay(atac_peakMatrix_pseudobulk.se[j,]))

###############################################################
## TEST: compare correlations before and after normalisation ##
###############################################################

# cor_before_norm.se <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/trajectories/nmp/CISBP_cor_TFexpr_vs_peakAcc.rds")
# cor_after_norm.se <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/results/rna_atac/rna_vs_acc/metacells/TFexpr_vs_peakAcc/trajectories/nmp/test/CISBP_cor_TFexpr_vs_peakAcc.rds")
# 
# cor_before_norm.cor <- as.numeric(assay(cor_before_norm.se,"cor"))
# cor_before_norm.cor <- cor_before_norm.cor[cor_before_norm.cor!=0]
# cor_after_norm.cor <- as.numeric(assay(cor_after_norm.se,"cor"))
# cor_after_norm.cor <- cor_after_norm.cor[cor_after_norm.cor!=0]
# 
# tmp <- sample.int(length(cor_before_norm.cor), size = 1e4)
# 
# to.plot <- rbind(
#   data.table(class = "before_norm", cor = cor_before_norm.cor[tmp]),
#   data.table(class = "after_norm", cor = cor_after_norm.cor[tmp])
# )
# 
# gghistogram(to.plot, x="cor", fill="class", bins=100) +
#   theme(
#     legend.title = element_blank(),
#     axis.text = element_text(size=rel(0.75), color="black")
#   )
# 
# atac_before_norm.mtx <- assay(atac_peakMatrix_metacells.se)
# # atac_after_norm.mtx <- log2(1e6*(sweep(atac_before_norm.mtx,2,colSums(atac_before_norm.mtx),"/"))+0.5)
# # atac_after_norm.mtx <- log2(1e6*(sweep(atac_before_norm.mtx,2,colSums(atac_before_norm.mtx),"/"))+0.5)
# atac_after_norm.mtx <- t(t(atac_before_norm.mtx) / colSums(atac_before_norm.mtx)) * 1e6
# atac_after_norm_log.mtx <- log2( (t(t(atac_before_norm.mtx) / colSums(atac_before_norm.mtx)) * 1e6)+1)
# 
# # i <- "T"; j <- "chr7:79789167-79789767"
# # i <- "GATA1"; j <- "chr3:137957674-137958274"
# # i <- "HOXA9"; j <- "chr18:42469870-42470470"
# i <- "FOXA2"; j <- "chr19:53135999-53136599"
# to.plot <- data.table(
#   gene = i,
#   peak = j,
#   cell = colnames(rna_metacells.sce),
#   celltype = metacell_metadata.dt$celltype,
#   expr = logcounts(rna_metacells.sce[i,])[1,],
#   raw = atac_before_norm.mtx[j,],
#   norm = atac_after_norm.mtx[j,],
#   norm_log = atac_after_norm_log.mtx[j,]
# ) %>% melt(id.vars=c("gene","peak","cell","celltype","expr"), variable.name="acc")
# 
# 
# ggplot(to.plot, aes(x=value, y=expr)) +
#   geom_point(aes(fill=celltype), shape=21) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   facet_wrap(~acc, scales = "free_x") +
#   labs(x="Chromatin accessibility (metacells)", y="Gene expression (metacells)") +
#   coord_cartesian(ylim=c(0,max(to.plot$expr)+0.25)) +
#   stat_smooth(method="lm", color="black") +
#   stat_cor(method = "pearson") +
#   theme_classic() +
#   theme(
#     axis.text = element_text(colour="black",size=rel(0.8)),
#     legend.position = "none"
#   )
