#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$archR.pseudobulk.GeneScoreMatrix.se <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn/pseudobulk_GeneScoreMatrix_TSS_summarized_experiment.rds")
io$archR.markers_genes <- file.path(io$basedir,"results_new/atac/archR/differential/GeneScoreMatrix_TSS/markers/marker_genes.txt.gz")
io$outdir <- file.path(io$basedir,"results_new/atac/archR/differential/GeneScoreMatrix_TSS/markers"); dir.create(io$outdir, showWarnings = F)

# Options
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  # "Mixed_mesoderm",
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
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

###################################
## Load precomputed marker genes ##
###################################

markers_genes.dt <- fread(io$archR.markers_genes) %>%
  .[celltype%in%opts$celltypes]

##########################
## Load pseudobulk ATAC ##
##########################

# Load SummarizedExperiment
atac.GeneScoreMatrix.se <- readRDS(io$archR.pseudobulk.GeneScoreMatrix.se)[,opts$celltypes]

# Create long data.table
atac_genes.dt <- assay(atac.GeneScoreMatrix.se) %>% as.data.table(keep.rownames = "name") %>%
  melt(id.vars="name", variable.name="celltype", value.name="acc")

################################################
## Plot number of marker genes per cell types ##
################################################

to.plot <- markers_genes.dt %>% .[,.N,by=c("celltype")]

p <- ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of marker genes") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none"
)

pdf(file.path(io$outdir, "barplot_number_marker_genes.pdf"), width = 9, height = 5)
print(p)
dev.off()

################################################
## Plot cell type exclusivity of marker genes ##
################################################

to.plot <- markers_genes.dt %>% .[,N:=.N,by="name"]

p <- ggboxplot(to.plot, x="celltype", y="N", fill="celltype", color="black", outlier.shape=NA) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Exclusivity of peak markers\n(the smaller the more exclusive)") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.title.y = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"boxplot_exclusivity_marker_genes.pdf"), width = 9, height = 5)
print(p)
dev.off()


################################################
## Plot cell type exclusivity of marker genes ##
################################################

to.plot <- markers_genes.dt %>%
  .[,.(Nx=.N),by="name"] %>%
  .[,Nx:=factor(Nx)] %>%
  .[,.(Ny=.N),by="Nx"]

p <- ggbarplot(to.plot, x="Nx", y="Ny", fill="gray70") +
  labs(x="Number of different cell types per marker peak", y="") +
  theme(
    axis.text = element_text(size=rel(0.75)),
  )

pdf(file.path(io$outdir,"boxplot_exclusivity2_marker_genes.pdf"), width = 7, height = 5)
print(p)
dev.off()



###################################################
## Plot chromatin accessibility por marker genes ##
###################################################

celltypes.to.plot <- c("Gut","Endothelium","Cardiomyocytes","Neural_crest")

# celltypes.to.plot <- opts$celltypes
opts$max.atac <- 1.5

to.plot <- atac_genes.dt %>%
  merge(markers_genes.dt[celltype%in%celltypes.to.plot,c("celltype","name")] %>% setnames("celltype","class"), by="name", allow.cartesian=T) %>%
  .[,class:=sprintf("%s markers",class)] %>%
  .[acc>=opts$max.atac,acc:=opts$max.atac] %>%
  # .[,.(expr=median(expr), acc=median(acc)), by=c("class","celltype")] %>%
  # .[,.(expr=mean(expr), acc=mean(acc)), by=c("class","celltype")] %>%
  .[,acc:=minmax.normalisation(acc)]

p <- ggplot(to.plot, aes_string(x="celltype", y="acc", fill="celltype")) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  facet_wrap(~class, ncol=2, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  coord_cartesian(ylim=c(0,0.60)) +
  labs(x="", y="Chromatin accessibility (scaled)") +
  theme_classic() +
  theme(
    axis.text.y = element_text(color="black", size=rel(0.75)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# pdf(file.path(io$outdir,"boxplot_peak_acc_celltype_distance_marker_genes.pdf"), width = 8, height = 6)
print(p)
# dev.off()

################################
## Distance from the epiblast ##
################################

sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]
reducedDim(sce.pseudobulk, "PCA") <- irlba::prcomp_irlba(t(logcounts(sce.pseudobulk)), n=10)$x

celltype_distance.dt <- dist(reducedDim(sce.pseudobulk, "PCA")) %>% as.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","celltypeA") %>%
  melt(id.vars="celltypeA", variable.name="celltypeB", value.name="distance")

to.plot <- atac_genes.dt %>%
  merge(markers_genes.dt[celltype=="Epiblast",c("celltype","name")] %>% setnames("celltype","class"), by="name", allow.cartesian=T) %>%
  .[,class:=sprintf("%s markers",class)] %>%
  .[acc>=opts$max.atac,acc:=opts$max.atac] %>%
  # .[,.(expr=median(expr), acc=median(acc)), by=c("class","celltype")] %>%
  # .[,.(expr=mean(expr), acc=mean(acc)), by=c("class","celltype")] %>%
  .[,acc:=minmax.normalisation(acc)] %>%
  .[,.(acc=mean(acc)), by=c("class","celltype")] %>%
  merge(
    celltype_distance.dt[celltypeA=="Epiblast"] %>% .[,celltypeA:=NULL] %>% setnames("celltypeB","celltype"), by="celltype"
  ) %>%
  .[,distance:=minmax.normalisation(distance)] 

tmp <- to.plot[celltype%in%c("Epiblast","Rostral_neurectoderm","Primitive_Streak","Somitic_mesoderm","NMP","Erythroid1")]


p <- ggscatter(to.plot, x="distance", y="acc", fill="celltype", size=4, shape=21,
          add="reg.line", add.params = list(color="gray40", fill="lightgray"), conf.int=TRUE, fullrange = TRUE) +
  stat_cor(method = "pearson") +
  # facet_wrap(~class, nrow=1, scales="free") +
  ggrepel::geom_text_repel(data=tmp, aes(x=distance, y=acc, label=celltype), size=4) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Transcriptomic distance from the Epiblast (scaled)", y="Chromatin accessibility of Epiblast markers") +
  theme(
    strip.text = element_text(size=rel(0.85)),
    axis.text = element_text(color="black", size=rel(0.70)),
    axis.title = element_text(color="black", size=rel(0.90)),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"marker_genes_epiblast_distance.pdf"), width = 8, height = 6)
print(p)
dev.off()
