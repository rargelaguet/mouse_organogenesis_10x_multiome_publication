# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$cell_metadata <- io$metadata
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/metacells_metadata.txt.gz")

io$cells_tf_expr_vs_gene_expr <- file.path(io$basedir,"results/rna/coexpression/correlation_matrix_tf2gene_single_cells.rds")
io$pseudobulk_tf_expr_vs_gene_expr <- file.path(io$basedir,"results/rna/coexpression/correlation_matrix_tf2gene_pseudobulk.rds")
io$metacells_tf_expr_vs_gene_expr <- file.path(io$basedir,"results/rna/coexpression/correlation_matrix_tf2gene_metacells.rds")

io$cells_sce <- io$rna.sce
io$metacells_sce <- file.path(io$basedir,"results/rna/metacells/SingleCellExperiment_metacells.rds")
io$pseudobulk_sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype.mapped//SingleCellExperiment_pseudobulk.rds")

io$outdir <- file.path(io$basedir,"results/rna/coexpression/pdf"); dir.create(io$outdir, showWarnings = F)
# Options

####################
## Load metadata  ##
####################

metacell_metadata.dt <- fread(io$metacell_metadata)

cell_metadata.dt <- fread(io$cell_metadata) %>% 
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & pass_atacQC==TRUE & doublet_call==FALSE & !is.na(celltype)]

########################################################
## Load TF expression vs gene expression correlations ##
########################################################

tf_expr_vs_gene_expr_cells.mtx <- readRDS(io$cells_tf_expr_vs_gene_expr)
tf_expr_vs_gene_expr_metacells.mtx <- readRDS(io$metacells_tf_expr_vs_gene_expr)
tf_expr_vs_gene_expr_pseudobulk.mtx <- readRDS(io$pseudobulk_tf_expr_vs_gene_expr)

opts$TFs <- rownames(tf_expr_vs_gene_expr_pseudobulk.mtx)

########################
## Load metacell data ##
########################

# Load RNA SingleCellExperiment
rna_metacells.sce <- readRDS(io$metacells_sce)
rna_metacells_tf.sce <- rna_metacells.sce[toupper(rownames(rna_metacells.sce))%in%opts$TFs]
rownames(rna_metacells_tf.sce) <- toupper(rownames(rna_metacells_tf.sce))

##########################
## Load pseudobulk data ##
##########################

# Load RNA SingleCellExperiment
rna_pseudobulk.sce <- readRDS(io$pseudobulk_sce)
rna_pseudobulk_tf.sce <- rna_pseudobulk.sce[toupper(rownames(rna_pseudobulk.sce))%in%opts$TFs]
rownames(rna_pseudobulk_tf.sce) <- toupper(rownames(rna_pseudobulk_tf.sce))

####################
## Load cell data ##
####################

# Load RNA SingleCellExperiment
rna_cells.sce <- load_SingleCellExperiment(file=io$cells_sce, cells=cell_metadata.dt$cell, remove_non_expressed_genes = T, normalise = T)
colData(rna_cells.sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame
rna_cells_tf.sce <- rna_cells.sce[toupper(rownames(rna_cells.sce))%in%opts$TFs]
rownames(rna_cells_tf.sce) <- toupper(rownames(rna_cells_tf.sce))

#########################
## Plot coverage stats ##
#########################

# Subset to marker genes
genes.to.plot <- unique(fread(io$rna.atlas.marker_genes.up)[["gene"]])
tfs.to.plot <- genes.to.plot[toupper(genes.to.plot)%in%rownames(tf_expr_vs_gene_expr_metacells.mtx)] %>% toupper
genes.to.plot <- genes.to.plot[genes.to.plot%in%colnames(tf_expr_vs_gene_expr_metacells.mtx)]

# tf_expr_vs_gene_expr_pseudobulk.cor <- tf_expr_vs_gene_expr_pseudobulk.mtx[tf_expr_vs_gene_expr_pseudobulk.mtx!=0]
# tf_expr_vs_gene_expr_metacells.cor <- tf_expr_vs_gene_expr_metacells.mtx[tf_expr_vs_gene_expr_metacells.mtx!=0]
# tf_expr_vs_gene_expr_cells.cor <- tf_expr_vs_gene_expr_cells.mtx[tf_expr_vs_gene_expr_cells.mtx!=0]
# tmp <- sample.int(length(tf_expr_vs_gene_expr_pseudobulk.mtx), size = 1e4)

to.plot <- rbind(
  data.table(class = "cell", cor = as.numeric(tf_expr_vs_gene_expr_cells.mtx[tfs.to.plot,genes.to.plot])),
  data.table(class = "metacell", cor = as.numeric(tf_expr_vs_gene_expr_metacells.mtx[tfs.to.plot,genes.to.plot])),
  data.table(class = "pseudobulk", cor = as.numeric(tf_expr_vs_gene_expr_pseudobulk.mtx[tfs.to.plot,genes.to.plot]))
)

p <- gghistogram(to.plot, x="cor", fill="gray70", bins=50) +
  facet_wrap(~class, nrow=1, scales="fixed") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(file.path(io$outdir,"tf_expr_vs_gene_expr_cor_histogram_cell_vs_metacell_vs_pseudobulk.pdf"), width=9, height=5)
print(p)
dev.off()

##########
## Plot ##
##########

TFs.to.plot <- c("GATA1")

# i <- "T"
for (i in TFs.to.plot) {
  
  tmp <- data.table(
    gene = genes.to.plot,
    pseudobulk_cor = tf_expr_vs_gene_expr_pseudobulk.mtx[i,genes.to.plot],
    metacells_cor = tf_expr_vs_gene_expr_metacells.mtx[i,genes.to.plot],
    cells_cor = tf_expr_vs_gene_expr_cells.mtx[i,genes.to.plot]
  ) %>% setorder(-pseudobulk_cor)
  
  genes.to.plot <- tmp$gene[-1] %>% head(n=5)
  
  # j <- genes.to.plot[1]
  for (j in genes.to.plot) {

    cells_expr.dt <- data.table(
      tf = i,
      gene = j,
      celltype = rna_cells.sce$celltype,
      tf_expr = logcounts(rna_cells_tf.sce[i,])[1,],
      gene_expr = logcounts(rna_cells.sce[j,])[1,]
    )

    p.cells <- ggplot(cells_expr.dt, aes(x=tf_expr, y=gene_expr)) +
      ggrastr::geom_point_rast(aes(color=celltype), size=0.5) +
      stat_smooth(method="lm", color="black") +
      scale_color_manual(values=opts$celltype.colors) +
      labs(x=sprintf("%s expression",i), y=sprintf("%s expression",j), title="Cells") +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(colour="black",size=rel(0.8)),
        legend.position = "none"
      )
    
    metacells_expr.dt <- data.table(
      tf = i,
      gene = j,
      celltype = rna_metacells.sce$celltype,
      tf_expr = logcounts(rna_metacells_tf.sce[i,])[1,],
      gene_expr = logcounts(rna_metacells.sce[j,])[1,]
    )
    
    p.metacells <- ggplot(metacells_expr.dt, aes(x=tf_expr, y=gene_expr)) +
      ggrastr::geom_point_rast(aes(color=celltype), size=1) +
      stat_smooth(method="lm", color="black") +
      scale_color_manual(values=opts$celltype.colors) +
      labs(x=sprintf("%s expression",i), y=sprintf("%s expression",j), title="Metacells") +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(colour="black",size=rel(0.8)),
        legend.position = "none"
      )
    
    pseudobulk_expr.dt <- data.table(
      tf = i,
      gene = j,
      celltype = colnames(rna_pseudobulk.sce),
      tf_expr = logcounts(rna_pseudobulk_tf.sce[i,])[1,],
      gene_expr = logcounts(rna_pseudobulk.sce[j,])[1,]
    )
    
    p.pseudobulk <- ggplot(pseudobulk_expr.dt, aes(x=tf_expr, y=gene_expr)) +
      geom_point(aes(fill=celltype), size=3, shape=21) +
      stat_smooth(method="lm", color="black") +
      scale_fill_manual(values=opts$celltype.colors) +
      labs(x=sprintf("%s expression",i), y=sprintf("%s expression",j), title="Pseudobulk") +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(colour="black",size=rel(0.8)),
        legend.position = "none"
      )
    
    p <- cowplot::plot_grid(plotlist=list(p.cells,p.metacells,p.pseudobulk), ncol=3)
    
    pdf(file.path(io$outdir,sprintf("%s_vs_%s_expr_cell_vs_metacell_vs_pseudobulk.pdf",i,j)), width=12, height=5.5)
    print(p)
    dev.off()
    
  }
}

