# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$basedir <- file.path(io$basedir,"test")
io$cell_metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
io$metacell_metadata <- file.path(io$basedir,"results/atac/archR/metacells/all_cells/PeakMatrix/metacells_metadata.txt.gz")
io$cell_atac_peak_matrix <- file.path(io$basedir,"processed/atac/archR/Matrices/PeakMatrix_summarized_experiment.rds")
io$metacell_atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/metacells/all_cells/PeakMatrix/PeakMatrix_summarized_experiment_metacells.rds")
io$pseudobulk_atac_peak_matrix <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype_genotype/PeakMatrix/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$outdir <- file.path(io$basedir,"results/atac/archR/plot_individual_peaks"); dir.create(io$outdir, showWarnings = F)

# Options
opts$samples <- c(
  "E8.5_CRISPR_T_KO",
  "E8.5_CRISPR_T_WT"
)

####################
## Load metadata  ##
####################

metacell_metadata.dt <- fread(io$metacell_metadata) %>% 
  .[sample%in%opts$samples & !is.na(celltype)]
table(metacell_metadata.dt$sample)

cell_metadata.dt <- fread(io$cell_metadata) %>% 
  .[pass_rnaQC==TRUE & pass_atacQC==TRUE & doublet_call==FALSE & !is.na(celltype) & sample%in%opts$samples]
table(cell_metadata.dt$sample)

##########################
## Load ATAC PeakMatrix ##
##########################

# atac_peak_matrix_cells.se <- readRDS(io$cell_atac_peak_matrix)[,cell_metadata.dt$cell]
atac_peak_matrix_metacells.se <- readRDS(io$metacell_atac_peak_matrix)[,metacell_metadata.dt$metacell]
atac_peak_matrix_pseudobulk.se <- readRDS(io$pseudobulk_atac_peak_matrix)

# Normalise ATAC data
# assay(atac_peak_matrix_metacells.se) <- t(t(assay(atac_peak_matrix_metacells.se)) / colSums(assay(atac_peak_matrix_metacells.se))) * 1e6
assay(atac_peak_matrix_metacells.se,"logcounts") <- log(1e6*(sweep(assay(atac_peak_matrix_metacells.se),2,colSums(assay(atac_peak_matrix_metacells.se),na.rm=T),"/"))+1)
assay(atac_peak_matrix_pseudobulk.se,"logcounts") <- log(1e6*(sweep(assay(atac_peak_matrix_pseudobulk.se),2,colSums(assay(atac_peak_matrix_pseudobulk.se),na.rm=T),"/"))+1)

# Sanity checks
stopifnot(rownames(atac_peak_matrix_metacells.se)==rownames(atac_peak_matrix_cells.se))
stopifnot(rownames(atac_peak_matrix_metacells.se)==rownames(atac_peak_matrix_pseudobulk.se))

###########
## Parse ##
###########

atac_peak_matrix_pseudobulk.se$celltype <- stringr::str_split(colnames(atac_peak_matrix_pseudobulk.se),"-") %>% map_chr(1)
atac_peak_matrix_pseudobulk.se$genotype <- stringr::str_split(colnames(atac_peak_matrix_pseudobulk.se),"-") %>% map_chr(2)

# Only consider cell types with sufficient number of metacells
stats.dt <- metacell_metadata.dt[,.N,by=c("celltype","genotype")] %>% dcast(celltype~genotype, value.var="N", fill=0)
celltypes.to.plot <- stats.dt[T_KO>=5 & WT>=5,celltype] 
stats.dt <- stats.dt[celltype%in%celltypes.to.plot]

# Subset data
# cell_metadata.dt <- cell_metadata.dt[celltype%in%celltypes.to.plot]
# atac_peak_matrix_cells.se <- atac_peak_matrix_cells.se[,cell_metadata.dt$cell]
metacell_metadata.dt <- metacell_metadata.dt[celltype%in%celltypes.to.plot]
atac_peak_matrix_metacells.se <- atac_peak_matrix_metacells.se[,metacell_metadata.dt$metacell]
atac_peak_matrix_pseudobulk.se <- atac_peak_matrix_pseudobulk.se[,atac_peak_matrix_pseudobulk.se$celltype%in%celltypes.to.plot]

##########
## Plot ##
##########

# stopifnot(colnames(atac_peak_matrix_cells.se)==cell_metadata.dt$cell)
stopifnot(colnames(atac_peak_matrix_metacells.se)==metacell_metadata.dt$metacell)

give.n <- function(x) { return(c(y = max(x), label = length(x)))}

features.to.plot <- sample(rownames(atac_peak_matrix_cells.se),10)
features.to.plot <- c("chr6:91880813-91881413")

# i <- features.to.plot[1]
for (i in features.to.plot) {
  
  ##################
  ## Prepare data ##
  ##################
  
  # single cells 
  cells_atac.dt <- data.table(
    cell = colnames(atac_peak_matrix_cells.se),
    # value = assay(atac_peak_matrix_cells.se)[i,],
    value = as.numeric(assay(atac_peak_matrix_cells.se)[i,]>0),
    celltype = cell_metadata.dt$celltype
  ) %>% .[,.(value=mean(value)), by="celltype"]# %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  # metacells
  metacells_atac.dt <- data.table(
    cell = colnames(atac_peak_matrix_metacells.se),
    value = assay(atac_peak_matrix_metacells.se,"logcounts")[i,],
    celltype = metacell_metadata.dt$celltype,
    genotype = metacell_metadata.dt$genotype
  )# %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  # pseudobulk
  pseudobulk_atac.dt <- data.table(
    value = assay(atac_peak_matrix_pseudobulk.se,"logcounts")[i,],
    celltype = atac_peak_matrix_pseudobulk.se$celltype,
    genotype = atac_peak_matrix_pseudobulk.se$genotype
  )# %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  # ylim <- max(c(max(pseudobulk_atac.dt$value), max(metacells_atac.dt$value), max(cells_atac.dt$value)))
  
  ##########
  ## Plot ##
  ##########
  
  p.cells <- ggplot(cells_atac.dt, aes(x=celltype, y=value, fill=celltype)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
    labs(x="", y="accessibility", title=i) +
    # coord_cartesian(ylim=c(0,ylim)) +
    # geom_text(data = tmp, aes(label=N), size=2.5) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      # axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  p.metacells <- ggplot(metacells_atac.dt, aes(x=celltype, y=value, fill=genotype)) +
    geom_violin(alpha=0.75) +
    geom_boxplot(outlier.shape=NA, alpha=0.75) +
    # stat_summary(fun.data = function(x) { return(c(y = max(metacells_atac.dt$value)+0.1, label = length(x)))}, geom = "text", size=2.75) +
    labs(x="", y="accessibility") +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(colour="black",size=rel(0.75)),
      # axis.text.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  p.pseudobulk <- ggplot(pseudobulk_atac.dt, aes(x=celltype, y=value, fill=genotype)) +
    geom_bar(stat="identity", color="black", position="dodge") +
    # scale_fill_manual(values=opts$celltype.colors[celltypes.to.plot]) +
    labs(x="", y="accessibility") +
    guides(x = guide_axis(angle = 90)) +
    # coord_cartesian(ylim=c(0,ylim)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(colour="black",size=rel(0.8)),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  # p <- cowplot::plot_grid(plotlist=list(p.cells,p.metacells,p.pseudobulk), ncol=1, rel_heights = c(0.3,0.3,0.4))
  p <- cowplot::plot_grid(plotlist=list(p.metacells,p.pseudobulk), ncol=2)
  
  pdf(file.path(io$outdir,sprintf("%s_cell_metacell_pseudobulk_acc.pdf",gsub("[:-]","_",i))), width=10, height=10)
  print(p)
  dev.off()
}


##########
## TEST ##
##########

i <- c("chr6:91880813-91881413")

# metacells
metacells_atac.dt <- data.table(
  cell = colnames(atac_peak_matrix_metacells.se),
  value = assay(atac_peak_matrix_metacells.se,"logcounts")[i,],
  celltype = metacell_metadata.dt$celltype,
  genotype = metacell_metadata.dt$genotype
)# %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]

# pseudobulk
pseudobulk_atac.dt <- data.table(
  value = assay(atac_peak_matrix_pseudobulk.se,"logcounts")[i,],
  celltype = atac_peak_matrix_pseudobulk.se$celltype,
  genotype = atac_peak_matrix_pseudobulk.se$genotype
)# %>% .[,celltype:=factor(celltype,levels=celltypes.to.plot)]

foo <- metacells_atac.dt[,mean(value),by=c("celltype","genotype")] %>% dcast(celltype~genotype, value.var="V1")
bar <- pseudobulk_atac.dt[,mean(value),by=c("celltype","genotype")] %>% dcast(celltype~genotype, value.var="V1")
