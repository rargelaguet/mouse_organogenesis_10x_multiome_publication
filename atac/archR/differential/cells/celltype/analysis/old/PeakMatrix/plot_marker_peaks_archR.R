#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$archR.peakSet.granges <- file.path(io$basedir,"processed_new/atac/archR/PeakCalls/PeakSet.rds")
io$archR.pseudobulk.peakMatrix.se <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$archR.markers_peaks <- file.path(io$basedir,"results_new/atac/archR/differential/PeakMatrix/markers/marker_peaks_v1.txt.gz")
io$outdir <- file.path(io$basedir,"results_new/atac/archR/differential/PeakMatrix/markers"); dir.create(io$outdir, showWarnings = F)

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

########################
## Load cell metadata ##
########################

# sample_metadata <- fread(io$metadata) %>%
#   .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
#   .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]
# 
# ########################
# ## Load ArchR project ##
# ########################
# 
# source(here::here("atac/archR/load_archR_project.R"))
# 
# ArchRProject.filt <- ArchRProject[sample_metadata$cell]
# 
# # add celltype.predicted to ArchR's CellColData
# tmp <- sample_metadata %>% 
#   .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
#   as.data.frame() %>% tibble::column_to_rownames("cell")
# stopifnot(all(tmp$TSSEnrichment_atac == getCellColData(ArchRProject.filt, "TSSEnrichment")[[1]]))
# ArchRProject.filt <- addCellColData(
#   ArchRProject.filt,
#   data = tmp$celltype.predicted, 
#   name = "celltype.predicted",
#   cells = rownames(tmp),
#   force = TRUE
# )

###################################
## Load precomputed marker genes ##
###################################

markers_peaks.dt <- fread(io$archR.markers_peaks) %>%
  .[celltype%in%opts$celltypes]

# Define GenomicRanges object with marker peaks
peakset.gr <- readRDS(io$archR.peakSet.granges)
peakset.gr$idx <- sprintf("%s:%s-%s",seqnames(peakset.gr), start(peakset.gr), end(peakset.gr))
peakset.gr <- peakset.gr[peakset.gr$idx %in% unique(markers_peaks.dt$idx)]

##########################
## Load pseudobulk ATAC ##
##########################

# Load SummarizedExperiment
atac.peakMatrix.se <- readRDS(io$archR.pseudobulk.peakMatrix.se)[,opts$celltypes]

# Load peak metadata
peak_metadata.dt <- fread(io$archR.peak.metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

# Define peak names
peak_names <- rowData(atac.peakMatrix.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac.peakMatrix.se) <- peak_names

# Subset peaks
peak_metadata.dt <- peak_metadata.dt[idx%in%unique(markers_peaks.dt$idx)]
atac.peakMatrix.se <- atac.peakMatrix.se[rownames(atac.peakMatrix.se) %in% peak_metadata.dt$idx]

# Create long data.table
atac_peaks.dt <- assay(atac.peakMatrix.se) %>% as.data.table(keep.rownames = "idx") %>%
  melt(id.vars="idx", variable.name="celltype", value.name="acc")

################################################
## Plot number of marker peaks per cell types ##
################################################

to.plot <- markers_peaks.dt %>% .[,.N,by=c("celltype")]

p <- ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of marker genes") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none"
)

pdf(file.path(io$outdir, "barplot_number_marker_peaks.pdf"), width = 9, height = 5)
print(p)
dev.off()

################################################
## Plot cell type exclusivity of marker peaks ##
################################################

to.plot <- markers_peaks.dt %>% .[,N:=.N,by="idx"]

p <- ggboxplot(to.plot, x="celltype", y="N", fill="celltype", color="black", outlier.shape=NA) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Exclusivity of peak markers\n(the smaller the more exclusive)") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.title.y = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"boxplot_exclusivity.pdf"), width = 9, height = 5)
print(p)
dev.off()


################################################
## Plot cell type exclusivity of marker peaks ##
################################################

to.plot <- markers_peaks.dt %>%
  .[,.(Nx=.N),by="idx"] %>%
  .[,Nx:=factor(Nx)] %>%
  .[,.(Ny=.N),by="Nx"]

p <- ggbarplot(to.plot, x="Nx", y="Ny", fill="gray70") +
  labs(x="Number of different cell types per marker peak", y="") +
  theme(
    axis.text = element_text(size=rel(0.75)),
  )

pdf(file.path(io$outdir,"boxplot_exclusivity2.pdf"), width = 7, height = 5)
print(p)
dev.off()



###################################################
## Plot chromatin accessibility por marker peaks ##
###################################################

celltypes.to.plot <- c("Gut","Endothelium","Cardiomyocytes","Neural_crest")

# celltypes.to.plot <- opts$celltypes
opts$max.atac <- 1.5

to.plot <- atac_peaks.dt %>%
  merge(markers_peaks.dt[celltype%in%celltypes.to.plot,c("celltype","idx")] %>% setnames("celltype","class"), by="idx", allow.cartesian=T) %>%
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

pdf(file.path(io$outdir,"boxplot_peak_acc_celltype_distance.pdf"), width = 8, height = 6)
print(p)
dev.off()


################################
## Distance from the epiblast ##
################################

sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]
reducedDim(sce.pseudobulk, "PCA") <- irlba::prcomp_irlba(t(logcounts(sce.pseudobulk)), n=10)$x

celltype_distance.dt <- dist(reducedDim(sce.pseudobulk, "PCA")) %>% as.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","celltypeA") %>%
  melt(id.vars="celltypeA", variable.name="celltypeB", value.name="distance")


to.plot <- atac_peaks.dt %>%
  merge(markers_peaks.dt[celltype=="Epiblast",c("celltype","idx")] %>% setnames("celltype","class"), by="idx", allow.cartesian=T) %>%
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

pdf(file.path(io$outdir,"epiblast_distance.pdf"), width = 8, height = 6)
print(p)
dev.off()
