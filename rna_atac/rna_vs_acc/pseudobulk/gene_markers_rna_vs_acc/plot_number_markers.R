#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$basedir <- file.path(io$basedir,"test")
io$marker_genes_rna <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_genes_filtered.txt.gz")
io$marker_peaks_atac <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype/PeakMatrix/parsed/markers_filt.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/gene_markers_rna_vs_acc")

# Options
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
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
  # "Mesenchyme",
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
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)


###############
## Load data ##
###############

markers_genes_rna.dt <- fread(io$marker_genes_rna) %>% .[celltype%in%opts$celltypes]
marker_peaks_atac.dt <- fread(io$marker_peaks_atac) %>% .[celltype%in%opts$celltypes]

##########
## Plot ##
##########

# Plot number of marker genes per cell types
to.plot <- rbind(
  markers_genes_rna.dt %>% .[,.N,by=c("celltype")] %>% .[,class:="Genes"],
  marker_peaks_atac.dt %>% .[,.N,by=c("celltype")] %>% .[,class:="ATAC peaks"]
) %>% .[,class:=factor(class, levels=c("Genes","ATAC peaks"))]

# Rename celltypes
opts$rename.celltypes <- c(
  "Forebrain_Midbrain_Hindbrain" = "Brain",
  "Haematoendothelial_progenitors" = "Haematoend. progenitors"
)
to.plot %>% .[,celltype:=stringr::str_replace_all(celltype,opts$rename.celltypes)] %>% .[,celltype:=gsub("_"," ",celltype)]
opts$celltype.colors["Haematoend. progenitors"] <- opts$celltype.colors["Haematoendothelial_progenitors"]
names(opts$celltype.colors) <- gsub("_"," ",names(opts$celltype.colors))

# Plot
p <- ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  facet_wrap(~class, ncol=2, scales="free_y") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of markers") +
  theme(
    strip.background = element_rect(colour="black", fill=NA),
    axis.text.y = element_text(size=rel(0.65)),
    axis.text.x = element_text(colour="black",size=rel(0.6), angle=90, hjust=1, vjust=0.5),
    axis.title = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_blank(),
    legend.position = "none"
)

pdf(file.path(io$outdir,"barplot_number_markers.pdf"), width = 7, height = 4)
print(p)
dev.off()

