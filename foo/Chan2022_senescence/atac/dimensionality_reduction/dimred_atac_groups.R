# here::i_am("dimensionality_reduction/dimensionality_reduction.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/atac/sample_metadata.txt")
# io$sce <- file.path(io$basedir,"results/atac/deeptools/SingleCellExperiment_chromHMM.rds")
# io$outdir <- file.path(io$basedir,"results/atac/dimensionality_reduction/groups/chromHMM"); dir.create(io$outdir, showWarnings = F)
io$sce <- file.path(io$basedir,"results/atac/deeptools/SingleCellExperiment_BJ_reprogramming.rds")
io$outdir <- file.path(io$basedir,"results/atac/dimensionality_reduction/groups/BJ_reprogramming"); dir.create(io$outdir, showWarnings = F)

# Options
opts$groups <- c(
  "p8_Control", 
  "p11_Control", 
  "p14_Control", 
  "p19_Control", 
  "p22_Control", 
  "p24_Control", 
  "p29_Control", 
  "p8_TG", 
  "p11_TG", 
  "p14_TG", 
  "p19_TG", 
  "p22_TG", 
  "p24_TG", 
  "p29_TG"
)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata) %>%
  .[Group%in%opts$groups,c("Group","Passage","Treatment")] %>% 
  .[,group:=sprintf("%s_p%s",Treatment,Passage)] %>%
  unique

###############
## Load data ##
###############

sce <- readRDS(io$sce)[,sample_metadata.dt$group]

# Add sample metadata as colData
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("group") %>% DataFrame

# Normalise
# assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(assay(sce,"counts")),"/"))+1)
assay(sce,"logcounts") <- assay(sce,"atac")

##################
## Subset genes ##
##################

# genes.to.use <- grep("^RP[L|S]",rownames(sce),value=T)
# genes.to.use <- grep("^MT-",rownames(sce),value=T)
# sce.filt <- sce[genes.to.use,]

#######################
## Feature selection ##
#######################

# decomp <- modelGeneVar(sce)
# hvgs <- rownames(decomp)[decomp$p.value<0.01 & decomp$mean > 1]
# sce.filt <- sce[hvgs,]
# dim(sce.filt)

tmp <- apply(assay(sce,"logcounts"),1,var)
sce.filt <- sce[tmp>0.01,]
sce.filt

#########
## PCA ##
#########

sce.filt <- runPCA(sce.filt, ncomponents=5)
reducedDims(sce) <- reducedDims(sce.filt)

# cor(sample_metadata.dt$age,reducedDim(sce.filt,"PCA"))

plotPCA(sce.filt, colour_by="Passage", shape_by="Treatment", ncomponents=c(1,2), point_alpha=1, point_size=5)
plotPCA(sce.filt, colour_by="10:107546169-107546689", ncomponents=c(1,2), point_alpha=1, point_size=5)
# plotPCA(sce.filt, colour_by="Passage", ncomponents=c(1,2), point_alpha=1, point_size=2.5)
# plotPCA(sce.filt, colour_by="Group", ncomponents=c(1,2), point_alpha=1, point_size=2.5)

var.explained <- attr(reducedDim(sce.filt,"PCA"),"varExplained")
var.explained <- (100*(var.explained/sum(var.explained))) %>% round(1)

##########
## Plot ##
##########

to.plot <- reducedDim(sce.filt) %>% as.data.table(keep.rownames = T) %>% 
  setnames("rn","group") %>%
  # melt(id.vars="sample", variable.name="pc") %>%
  merge(sample_metadata.dt,by="group")

p <- ggplot(to.plot, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill=Passage, shape=Treatment), shape=21, size=5, stroke=0.35) +
  # ggrepel::geom_text_repel(aes_string(label="surface_marker"), data=to.plot[,.(PC1=median(PC1), PC2=median(PC2)), by="surface_marker"]) +
  scale_fill_viridis_c() +
  theme_classic() +
  labs(x=sprintf("PC1 (%.2f%%)",var.explained[1]), y=sprintf("PC1 (%.2f%%)",var.explained[2])) +
  # ggplot_theme_NoAxes() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color="black")
    # axis.ticks = element_blank()
  )

pdf(file.path(io$outdir,"pca_Chan2022_senescence_atac.pdf"), width=8, height=6)
print(p)
dev.off()
