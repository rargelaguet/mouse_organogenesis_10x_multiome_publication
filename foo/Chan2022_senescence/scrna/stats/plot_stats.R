# here::i_am("rna/dimensionality_reduction/dimensionality_reduction_sce.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"cell_metadata.txt.gz")
io$sce <- file.path(io$basedir,"processed/cellranger/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results/sample_stats"); dir.create(io$outdir, showWarnings = F)

# Options
opts$samples <- c( 
  # "hTert_t2",
  # "hTert_t4",
  # "hTert_t5",
  # "hTert_t6",
  # "hTert_t7",
  "PDL_25_rep1",
  "PDL_29_rep1",
  "PDL_29_rep2",
  "PDL_33_rep1",
  "PDL_33_rep2",
  "PDL_37_rep1",
  "PDL_37_rep2",
  "PDL_46_rep1",
  "PDL_46_rep2",
  "PDL_50_rep1",
  "PDL_50_rep2"
)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[,sample:=NULL] %>% setnames("alias","sample") %>%
  .[sample%in%opts$samples] %>% .[,sample:=factor(sample,levels=opts$samples)]# %>%
# .[,stage:=factor(stage, levels=opts$stages, ordered = T)]

cell_metadata.dt[grep("PDL_50",sample),c("cell_cycle_s_score","cell_cycle_g2m_score","cell_cycle_phase"):=list(NA,NA,"G1")]

table(cell_metadata.dt$sample)
table(cell_metadata.dt$genetic_background)

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame


############################
## Plot cell cycle states ##
############################

to.plot <- cell_metadata.dt[,N:=.N,by=c("sample","group")] %>%
  .[,.(prop=100*.N/unique(N)),by=c("sample","group","cell_cycle_phase")]

p <- ggplot(to.plot, aes(x=sample, y=prop)) +
  geom_bar(aes(fill=cell_cycle_phase), stat="identity", color="black") +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  labs(x="", y="Cell cycle phase (%)") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    strip.text = element_text(color="black", size=rel(1.05)),
    axis.text.y = element_text(color="black", size=rel(1.15)),
    axis.text.x = element_text(color="black", size=rel(1)),
    axis.title = element_text(color="black", size=rel(1.0))
  )

pdf(file.path(io$outdir,"cell_cycle_proportions_Chan2022_senescence_scrna.pdf"), width=9, height=5)
print(p)
dev.off()