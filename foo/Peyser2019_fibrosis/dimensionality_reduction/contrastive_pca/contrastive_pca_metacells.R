# here::i_am("dimensionality_reduction/dimensionality_reduction_sce.R")

library(scran)
library(scater)

source(here::here("settings.R"))
source(here::here("utils.R"))

set.seed(42)

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/metacells/per_condition/metacells_metadata.txt.gz")
io$sce <- file.path(io$basedir,"results/metacells/per_condition/SingleCellExperiment_metacells.rds")
io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction/contrastive_PCA"); dir.create(io$outdir, showWarnings = F)

# Options
# opts$conditions <- c("Saline_Control", "Bleomycin_day11", "Bleomycin_nintedanib_day11")
opts$conditions <- c("Saline_Control", "Bleomycin_day11")

opts$celltypes <- c("dendritic_cells", "Type_II_alveolar", "Type_I_alveolar", "mesothelium", 
                    "fibroblasts", "macrophages", "endothelial", "B_cells", "monocytes", "neutrophils")

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  setnames("metacell","cell") %>%
  .[condition%in%opts$conditions] %>% .[,condition:=factor(condition,levels=opts$conditions)] %>%
  .[cell_type%in%opts$celltypes] %>% .[,cell_type:=factor(cell_type,levels=opts$celltypes)]

table(cell_metadata.dt$condition)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=2500) %>% rownames
sce_filt <- sce[hvgs,]

#####################
## contrastive PCA ##
#####################

# groups.target <- c("D00_Dox_Recover", "D02_Dox_Recover", "D04_Dox_Recover", "D06_Dox_Recover", "D08_Dox_Recover")
# groups.background <- c("D00_noDox", "D04_noDox", "D08_noDox", "D12_noDox")

groups.target <- "Bleomycin_day11"
groups.background <- "Saline_Control"

cPCA.out <- scPCA(
  target = t(assay(sce[,sce$group%in%groups.target], "logcounts")),
  background = t(assay(sce[,sce$group%in%groups.background], "logcounts")),
  alg = "var_proj",
  penalties = 0, # no penalties = non-sparse cPCA.
  n_eigen = 15,
  contrasts = 10
)

reducedDim(target, "cPCA") <- cPCA.out$x


##########
## Plot ##
##########

dimred.to.plot <- "PCA"
# dimred.to.plot <- "cPCA"

to.plot <- reducedDim(sce,dimred.to.plot)[,c(1,2)] %>% as.data.table %>% 
  setnames(c("V1","V2")) %>%
  .[,cell:=colnames(sce)] %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(fill=group), shape=21, size=1.50, stroke=0.05, color="black", alpha=0.90, raster.dpi=300) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  # facet_wrap(~treatment, nrow=2) +
  # ggrepel::geom_text_repel(aes(label=celltype), data=to.plot[,.(V1=median(V1), V2=median(V2)), by=c("celltype")]) +
  # ggrepel::geom_text_repel(aes(label=celltype), data=to.plot[,.(V1=median(V1), V2=median(V2)), by=c("celltype","treatment")]) +
  # scale_fill_manual(values=opts$colors) +
  labs(x=sprintf("%s-1",dimred.to.plot), y=sprintf("%s-2",dimred.to.plot)) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.9), color="black"),    
    legend.position  ="right",
    legend.title=element_blank()
  )

# pdf(file.path(io$outdir,"BJ_CHAMP_UMAP_celltype_metacells.pdf"), width=6, height=4.5)
# print(p)
# dev.off()

