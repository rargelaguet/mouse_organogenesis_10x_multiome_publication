here::i_am("rna/plot_individual_genes/genotype/plot_individual_genes_by_genotype.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$basedir <- file.path(io$basedir,"test")
io$cell_metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/all_cells/metacells_metadata.txt.gz")
io$rna_cells.sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
io$rna_metacells.sce <- file.path(io$basedir,"results/rna/metacells/all_cells/SingleCellExperiment_metacells.rds")
io$rna_pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype_genotype/SingleCellExperiment_pseudobulk.rds")
io$outdir <- file.path(io$basedir,"results/rna/plot_individual_genes/genotype"); dir.create(io$outdir, showWarnings = F, recursive = T)

opts$samples <- c(
  "E8.5_CRISPR_T_KO",
  "E8.5_CRISPR_T_WT"
)

opts$celltypes = c(
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

opts$min.cells <- 25

###################
## Load metadata ##
###################

# Cells
cell_metadata.dt <- fread(io$cell_metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples & celltype%in%opts$celltypes] %>%
  # setnames("celltype.mapped","celltype") %>%
  .[,celltype:=factor(celltype, levels=opts$celltypes)]

# Metacells
metacell_metadata.dt <- fread(io$metacell_metadata) %>% 
  .[sample%in%opts$samples & celltype%in%opts$celltypes] %>%
  .[,celltype:=factor(celltype, levels=opts$celltypes)]


#########################
## Load RNA expression ##
#########################

# Cells
sce <- load_SingleCellExperiment(io$rna_cells.sce, cells=cell_metadata.dt$cell, normalise = TRUE)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

# Metacells
metacells.sce <- readRDS(io$rna_metacells.sce)[,metacell_metadata.dt$metacell]
colData(metacells.sce) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

# Pseudobulk
pseudobulk.sce <- load_SingleCellExperiment(io$rna_pseudobulk.sce)
pseudobulk.sce$celltype <- names(pseudobulk.sce@metadata$n_cells) %>% strsplit(split="-") %>% map_chr(1)
pseudobulk.sce$genotype <- names(pseudobulk.sce@metadata$n_cells) %>% strsplit(split="-") %>% map_chr(2)

##########
## Plot ##
##########

# genes.to.plot <- c("Cyp26a1","Rspo3","Fgf8","Fgf17","Wnt3a","Limch1","Hes7","Mgst1")
genes.to.plot <- c("")
# genes.to.plot <- fread(io$atlas.marker_genes) %>% .[grep("Erythroid",celltype),gene] %>% unique

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]
celltypes.to.plot <- c("Spinal_cord","NMP","Caudal_Mesoderm","Somitic_mesoderm")
celltypes.to.plot <- opts$celltypes

i <- "Hoxc9"
# Cped1
for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce)[i,],
    genotype = sce$genotype,
    celltype = sce$celltype
  ) %>% .[celltype%in%celltypes.to.plot] %>%
    .[,N:=.N,by=c("genotype","celltype")] %>% .[N>=opts$min.cells] %>%
    .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  p <- ggplot(to.plot, aes(x=genotype, y=expr, fill=genotype)) +
    geom_violin(scale = "width", alpha=0.8) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
    stat_summary(fun.data = give.n, geom = "text", size=3) +
    # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
    # scale_fill_manual(values=opts$colors) +
    scale_fill_brewer(palette="Dark2") +
    facet_wrap(~celltype, scales="fixed") +
    theme_classic() +
    labs(x="",y=sprintf("%s expression",i)) +
    theme(
      strip.text = element_text(size=rel(0.85)),
      # axis.text.x = element_text(colour="black",size=rel(0.9)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(0.85))
    )
  
  pdf(sprintf("%s/%s_boxplots_single_cells_wt_vs_ko.pdf",io$outdir,i), width=10, height=9)
  print(p)
  dev.off()

}

####################
## Plot metacells ##
####################

# genes.to.plot <- c("Cyp26a1","Rspo3","Fgf8","Fgf17","Wnt3a","Limch1","Hes7","Mgst1")
genes.to.plot <- c("T")

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]
celltypes.to.plot <- c("Spinal_cord","NMP","Caudal_Mesoderm","Somitic_mesoderm")

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    cell = colnames(metacells.sce),
    expr = logcounts(metacells.sce)[i,],
    genotype = metacells.sce$genotype,
    celltype = metacells.sce$celltype
  ) %>% .[celltype%in%celltypes.to.plot] %>%
    .[,N:=.N,by=c("genotype","celltype")] %>%
    .[,celltype:=factor(celltype,levels=celltypes.to.plot)]
  
  p <- ggplot(to.plot, aes(x=genotype, y=expr, fill=genotype)) +
    geom_violin(scale = "width", alpha=0.8) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
    stat_summary(fun.data = give.n, geom = "text", size=3) +
    # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
    # scale_fill_manual(values=opts$colors) +
    scale_fill_brewer(palette="Dark2") +
    facet_wrap(~celltype, scales="fixed") +
    theme_classic() +
    labs(x="",y=sprintf("%s expression",i)) +
    theme(
      strip.text = element_text(size=rel(0.85)),
      # axis.text.x = element_text(colour="black",size=rel(0.9)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(0.85))
    )
  
  pdf(sprintf("%s/%s_boxplots_metacells_wt_vs_ko.pdf",io$outdir,i), width=10, height=9)
  print(p)
  dev.off()
  
}

################################
## Plot heatmaps (pseudobulk) ##
################################

# genes.to.plot <- "T"
# 
# genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]
# celltypes.to.plot <- c("Spinal_cord","NMP","Caudal_mesoderm","Somitic_mesoerm")
# 
# for (i in i) {
#   
#   tmp <- data.table(
#     expr = logcounts(pseudobulk.sce)[gene,],
#     genotype = pseudobulk.sce$genotype,
#     celltype = pseudobulk.sce$celltype
#   ) %>% .[celltype%in%celltypes.to.plot]
#   
#   to.plot <- data.table(expand.grid(opts$celltypes,unique(pseudobulk.sce$genotype))) %>%
#     setnames(c("celltype","genotype")) %>% 
#     merge(cell_metadata.dt[,.N,c("celltype","genotype")],all.x=TRUE) %>%
#     merge(tmp,by=c("celltype","genotype"),all.x=TRUE) %>%
#     .[N<=30,expr:=NA]
# 
#   p <- ggplot(to.plot, aes(x=celltype, y=genotype, fill=expr)) +
#     geom_tile(color="black") +
#     # geom_text(aes(label=round(expr,1))) +
#     scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = 'gray70') +
#     theme_classic() +
#     guides(x = guide_axis(angle = 90)) +
#     theme(
#       axis.text = element_text(color="black", size=rel(0.7)),
#       axis.title = element_blank(),
#       strip.background = element_blank(),
#       axis.ticks = element_blank(),
#       axis.line = element_blank(),
#       legend.title = element_blank()
#     )
#   
#   pdf(sprintf("%s/%s_heatmap_pseudobulk_wt_vs_ko.pdf",io$outdir,i), width=10, height=4)
#   print(p)
#   dev.off()
#   
# }

################################
## Plot barplots (pseudobulk) ##
################################

genes.to.plot <- "T"
# genes.to.plot <- c("Cyp26a1","Rspo3","Fgf8","Fgf17","Wnt3a","Limch1","Hes7","Mgst1")

genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]
# celltypes.to.plot <- cell_metadata.dt[,.N,by=c("celltype","genotype")] %>% .[N>=50] %>% .[,.N,by="celltype"] %>% .[N==2,celltype] %>% unique %>% as.character
celltypes.to.plot <- c("Spinal_cord","NMP","Caudal_mesoderm","Somitic_mesoerm")
celltypes.to.plot <- opts$celltypes

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    expr = logcounts(pseudobulk.sce)[i,],
    genotype = pseudobulk.sce$genotype,
    celltype = pseudobulk.sce$celltype
  ) %>% .[celltype%in%celltypes.to.plot]
  
  p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=genotype)) +
    geom_bar(stat="identity", position=position_dodge(width=0.7), color="black", width=0.7) +
    theme_classic() +
    guides(x = guide_axis(angle = 90)) +
    labs(x="", y="RNA expression") +
    theme_classic() +
    theme(
      axis.text.x = element_text(color="black", size=rel(0.75)),
      axis.text.y = element_text(color="black", size=rel(1))
    )
  
  pdf(sprintf("%s/%s_barplots_pseudobulk_wt_vs_ko.pdf",io$outdir,i), width=8, height=5)
  print(p)
  dev.off()
  
}

##########
## TEST ##
##########

# i <- "Epha5"
# j <- "ExE_endoderm"

# tmp <- sce[i,sce$celltype==j]
# mean(logcounts(tmp)[1,][tmp$genotype=="WT"])
# mean(logcounts(tmp)[1,][tmp$genotype=="T_KO"])

# tmp <- pseudobulk.sce[i,pseudobulk.sce$celltype==j]
# mean(logcounts(tmp)[1,][tmp$genotype=="WT"])
# mean(logcounts(tmp)[1,][tmp$genotype=="T_KO"])
# logcounts(tmp[i,])

# tmp <- metacells.sce[i,metacells.sce$celltype==j]
# # counts(tmp)[1,]
# mean(logcounts(tmp)[1,][tmp$genotype=="WT"])
# mean(logcounts(tmp)[1,][tmp$genotype=="T_KO"])

