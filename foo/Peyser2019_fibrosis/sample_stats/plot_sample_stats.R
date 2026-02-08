# here::i_am("rna/dimensionality_reduction/dimensionality_reduction_sce.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"processed/cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/sample_stats"); dir.create(io$outdir, showWarnings = F)

# Options
# opts$conditions <- c("Saline Control", "Bleomycin Fibrosis Model", "Bleomycin Fibrosis plus Nintedanib")

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] %>% .[,condition:=factor(condition,levels=opts$conditions)] %>%
  .[cell_type%in%opts$celltypes] %>% .[,cell_type:=factor(cell_type,levels=opts$celltypes)]

###############
## Load data ##
###############

# # Load RNA expression data as SingleCellExperiment object
# sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = FALSE)
# 
# # Add sample metadata as colData
# colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

############################
## Plot cell cycle phases ##
############################

# to.plot <- cell_metadata.dt[,.N,by=c("sample","cell_cycle_phase")] %>%
#   .[,fraction_cells:=N/sum(N),by="sample"]
# 
# p <- ggplot(to.plot, aes(x=sample, y=fraction_cells, fill=cell_cycle_phase)) +
#   geom_bar(stat="identity", color="black") +
#   labs(x="", y="Fraction of cells") +
#   theme_classic() +
#   theme(
#     # axis.text = element_text(color="black", size=rel(1.15)),
#     legend.title = element_blank(),
#     legend.position = "top"
#   )
# 
# pdf(file.path(io$outdir,sprintf("cell_cycle_phases_fibros.pdf")), width=7, height=6)
# print(p)
# dev.off()

#############################################
## Plot celltype proportions per condition ##
#############################################

to.plot <- cell_metadata.dt %>% .[,ncells:=.N,by=c("condition","sample")] %>% .[,.(N=.N, fraction_cells=.N/unique(ncells)),by=c("cell_type","condition","sample")]
# to.plot <- cell_metadata.dt[,.N,by=c("Tissue","ClustName")]

p <- ggplot(to.plot, aes(x=condition, y=N, fill=condition)) +
  # geom_bar(stat="identity", color="black") +
  geom_point(color="black", size=2.5, shape=21) +
  facet_wrap(~cell_type, scales="free_y") +
  # facet_grid(~condition, scales="free_x", space="free_x") +
  labs(x="", y="Number of cells") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size=rel(1.0), color="black"),
    # axis.text.x = element_text(colour="black",size=rel(0.85), angle=90, hjust=1, vjust=0.5),
    axis.text.x = element_blank(),
    # axis.title = element_text(colour="black",size=rel(0.75)),
    legend.position = "top",
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,sprintf("ncells_per_condition_Peyser2019_fibrosis.pdf")), width=8, height=7)
print(p)
dev.off()
