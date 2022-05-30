
# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$housekeeping.genes <-"/Users/argelagr/data/genesets/manual_genesets/housekeeping/housekeeping.tsv"
# io$marker_genes_rna <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_genes_filtered.txt.gz")
# io$rna.pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
# io$archR.pseudobulk.GeneMatrix.se <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype/GeneScoreMatrix_TSS/pseudobulk_GeneScoreMatrix_TSS_summarized_experiment.rds")
io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/peak_markers_rna_vs_acc"); dir.create(io$outdir, showWarnings = F)

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
  # "ExE_ectoderm"
  # "Parietal_endoderm"
)

####################
## Load gene sets ##
####################

# Housekeeping genes
houseekeping.dt <- fread(io$housekeeping.genes, header=F) %>% setnames(c("ens_id","gene"))

# Marker peaks
marker_peaks.dt <- fread(io$atac.markers_peaks.pseudobulk.filt) %>% .[celltype%in%opts$celltypes]

# Marker genes
marker_genes.dt <- fread(io$rna.celltype.marker_genes.filt) %>% .[celltype%in%opts$celltypes]

####################
## Load peak metadata ##
####################

peak_metadata.dt <- fread(io$archR.peak.metadata) %>% 
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[peak%in%unique(marker_peaks.dt$feature)]

# Load links between peaks and genes
peak2gene.dt <- fread(io$archR.peak2gene.nearest) %>% .[peak%in%unique(marker_peaks.dt$feature)]

##########################
## Load pseudobulk ATAC ##
##########################

# Load SummarizedExperiment
atac_GeneScores.se <- readRDS(io$archR.GeneScoreMatrix.pseudobulk)[,opts$celltypes]
atac_peakMatrix.se <- readRDS(io$archR.peakMatrix.pseudobulk)[unique(marker_peaks.dt$feature),opts$celltypes]

# Subset to non-promoter peaks
atac_peakMatrix.se <- atac_peakMatrix.se[rownames(atac_peakMatrix.se)%in%peak_metadata.dt[peakType!="Promoter",peak]]

# Normalise ATAC data
assay(atac_peakMatrix.se,"logcounts") <- log(1e6*(sweep(assay(atac_peakMatrix.se),2,colSums(assay(atac_peakMatrix.se),na.rm=T),"/"))+1)
assay(atac_GeneScores.se,"logcounts") <- log(1e6*(sweep(assay(atac_GeneScores.se),2,colSums(assay(atac_GeneScores.se),na.rm=T),"/"))+1)

# Prepare data.table
atac_genes.dt <- assay(atac_GeneScores.se,"logcounts") %>% as.data.table(keep.rownames = "gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="acc")
atac_peaks.dt <- assay(atac_peakMatrix.se,"logcounts") %>% as.data.table(keep.rownames = "feature") %>%
  melt(id.vars="feature", variable.name="celltype", value.name="acc")

# Cap outlier values
opts$min.atac <- 1
opts$max.atac <- 5.5
atac_genes.dt %>% .[acc<=opts$min.atac,acc:=opts$min.atac] %>% .[acc>=opts$max.atac,acc:=opts$max.atac]
atac_peaks.dt %>% .[acc<=opts$min.atac,acc:=opts$min.atac] %>% .[acc>=opts$max.atac,acc:=opts$max.atac]
  
#########################
## Load pseudobulk RNA ##
#########################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

# Prepare data.table
rna.dt <- logcounts(rna_pseudobulk.sce) %>%
  as.data.table(keep.rownames="gene") %>% 
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

# Cap outlier values
opts$min.expr <- 0.5
opts$max.expr <- 12
rna.dt %>% .[expr<=opts$min.expr,expr:=opts$min.expr] %>% .[expr>=opts$max.expr,expr:=opts$max.expr]
  
# Split gene sets 
rna_housekeeping.dt <- rna.dt[gene%in%houseekeping.dt$gene] %>%
  .[,class:="Positive control (housekeeping genes)"]

rna_olfactory.dt <- rna.dt[grep("Olfr",gene)] %>% 
  .[,class:="Negative control (olfactory receptors)"]

rna_markers.dt <- rna.dt %>%
  .[gene%in%marker_genes.dt[,gene]] %>%
  .[,class:="Marker genes"]

######################################################
## Scatterplot of peak vs gene accessibility levels ##
######################################################

celltypes.to.plot <- opts$celltypes
markers.to.plot <- c("Gut","Endothelium","Cardiomyocytes","Neural_crest")

atac_peaks_tmp.dt <- atac_peaks.dt %>% 
  .[celltype%in%celltypes.to.plot] %>%
  merge(marker_peaks.dt[celltype%in%markers.to.plot,c("celltype","feature")] %>% setnames("celltype","class"), by="feature", allow.cartesian=TRUE) %>%
  .[,.(mean=mean(acc), sd=sd(acc)),by=c("celltype","class")] %>%
  .[,feature:="Peaks"]

atac_genes_tmp.dt <- atac_genes.dt %>% 
  .[celltype%in%celltypes.to.plot] %>%
  setnames("gene","feature") %>%
  merge(marker_genes.dt[celltype%in%markers.to.plot,c("celltype","gene")] %>% setnames(c("class","feature")), by="feature", allow.cartesian=TRUE) %>%
  .[,.(mean=mean(acc), sd=sd(acc)),by=c("celltype","class")] %>%
  .[,feature:="Genes"]
  
to.plot <- rbindlist(list(atac_peaks_tmp.dt,atac_genes_tmp.dt))

p <- ggplot(to.plot, aes_string(x="celltype", y="mean", shape="feature", fill="celltype")) +
  geom_point(size=3.5) + 
  geom_line(aes(group=celltype, color=celltype)) +
  facet_wrap(~class, nrow=2, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  scale_color_manual(values=opts$celltype.colors) +
  scale_shape_manual(values=c(21,22))  +
  # coord_cartesian(ylim=c(0,opts$max.atac)) +
  labs(x="", y="Chromatin accessibility levels") +
  theme_classic() +
  theme(
    axis.text.y = element_text(color="black", size=rel(0.9)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"marker_genes_accessibility_scatterplots_v1.pdf"), width = 7, height = 5)
print(p)
dev.off()


to.plot2 <- to.plot %>% dcast(celltype+class~feature,value.var="mean")

p <- ggplot(to.plot2, aes_string(x="Genes", y="Peaks", fill="celltype")) +
  geom_point(shape=21, size=4) + 
  geom_abline(slope=1, intercept=0) +
  facet_wrap(~class, nrow=2, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  # scale_color_manual(values=opts$celltype.colors) +
  coord_cartesian(xlim = c(1,4), ylim = c(1,4)) +
  labs(x="Chromatin accessibility (Genes)", y="Chromatin accessibility (Peaks)") +
  theme_classic() +
  theme(
    axis.text = element_text(color="black", size=rel(0.9)),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"marker_genes_accessibility_scatterplots_v2.pdf"), width = 6.5, height = 5)
print(p)
dev.off()

##################################################################################
## Scatterplot of average RNA expression versus average chromatin accessibility ##
##################################################################################

celltypes.to.plot <- c("Gut","Endothelium","Cardiomyocytes","Neural_crest")
celltypes.to.plot <- opts$celltypes

for (i in celltypes.to.plot) {
  
  tmp <- peak2gene.dt[peak%in%marker_peaks.dt[celltype==i,feature] & gene%in%marker_genes.dt[celltype==i,gene]]
  peaks.to.use <- unique(tmp$peak)
  genes.to.use <- unique(tmp$gene)
  
  rna_houseekeping_avg.dt <- rna_housekeeping.dt %>% .[,.(expr=mean(expr)), by=c("celltype")] %>% .[,class:="Housekeeping genes"] 
  rna_olfactory_avg.dt <- rna_olfactory.dt %>% .[,.(expr=mean(expr)), by=c("celltype")] %>% .[,class:="Olfactory receptors"] 
  rna_markers_avg.dt <- rna_markers.dt[gene%in%genes.to.use] %>% .[,.(expr=mean(expr)), by=c("celltype")] %>% .[,class:="Marker genes"] 
  
  atac_houseekeping_genes_avg.dt <- atac_genes.dt[gene%in%houseekeping.dt$gene] %>% .[,.(acc=mean(acc)), by=c("celltype")] %>% .[,class:="Housekeeping genes"]
  atac_olfactory_genes_avg.dt <- atac_genes.dt[grep("Olfr",gene)] %>% .[,.(acc=mean(acc)), by=c("celltype")] %>% .[,class:="Olfactory receptors"]
  atac_marker_peaks_avg.dt <- atac_peaks.dt[feature%in%peaks.to.use] %>% .[,.(acc=mean(acc)), by=c("celltype")] %>% .[,class:="Marker peaks"]
  atac_marker_genes_avg.dt <- atac_genes.dt[gene%in%genes.to.use] %>% .[,.(acc=mean(acc)), by=c("celltype")] %>% .[,class:="Marker genes"]
  
  rna_atac_olfactory_genes_avg.dt <- merge(rna_olfactory_avg.dt, atac_olfactory_genes_avg.dt, by=c("celltype","class"))
  rna_atac_housekeeping_genes_avg.dt <- merge(rna_houseekeping_avg.dt, atac_houseekeping_genes_avg.dt, by=c("celltype","class"))
  rna_atac_marker_genes_avg.dt <- merge(rna_markers_avg.dt, atac_marker_genes_avg.dt, by=c("celltype","class"))
  rna_atac_marker_peaks_avg.dt <- merge(rna_markers_avg.dt %>% copy %>% .[,class:="Marker peaks"] , atac_marker_peaks_avg.dt, by=c("celltype","class"))
  
  to.plot <- rbind(rna_atac_olfactory_genes_avg.dt,rna_atac_housekeeping_genes_avg.dt,rna_atac_marker_genes_avg.dt,rna_atac_marker_peaks_avg.dt) %>%
    .[acc<=opts$min.atac,acc:=opts$min.atac] %>% .[expr<=opts$min.expr,expr:=opts$min.expr] %>%
    .[acc>=opts$max.atac,acc:=opts$max.atac] %>% .[expr>=opts$max.expr,expr:=opts$max.expr] %>%
    .[,c("expr","acc"):=list(minmax.normalisation(expr),minmax.normalisation(acc))]
    
  p <- ggplot(to.plot, aes_string(x="acc", y="expr")) +
    geom_jitter(aes(fill=celltype), shape=21, size=3, width=0.03, height=0.03) +
    geom_hline(yintercept=0.5, linetype="dashed") +
    geom_vline(xintercept=0.5, linetype="dashed") +
    scale_shape_manual(values=c(22,24)) +
    scale_fill_manual(values=opts$celltype.colors) + 
    labs(x="Chromatin accessibility (scaled)", y="RNA expression (scaled)") +
    coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
    scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
    facet_wrap(~class) +
    theme_classic() +
    theme(
      axis.text = element_text(color="black", size=rel(0.75)),
      axis.title = element_text(color="black", size=rel(1)),
      legend.position = "none",
      legend.title = element_blank()
    )
  
  pdf(file.path(io$outdir,sprintf("%s_marker_peaks_rna_vs_acc.pdf",i)), width = 7, height = 4)
  print(p)
  dev.off()
}

