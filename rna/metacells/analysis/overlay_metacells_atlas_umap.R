# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

source(here::here("rna/mapping/analysis/plot_utils.R"))

#####################
## Define settings ##
#####################

## I/O
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/metacells_metadata.txt.gz")
io$metacell_sce <- file.path(io$basedir,"results/rna/metacells/SingleCellExperiment_metacells.rds")
# io$umap <- file.path(io$basedir,"results/rna/dimensionality_reduction/sce/batch_correction_by_sample_remove_ExE_cells_False/umap_features2500_pcs50_neigh25_dist0.5.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/metacells/pdf"); dir.create(io$outdir, showWarnings = F)

# Dot size
opts$size.mapped <- 0.30
opts$size.nomapped <- 0.1

# Dot transparency
opts$alpha.mapped <- 0.75
opts$alpha.nomapped <- 0.35

###################
## Load metadata ##
###################

metacell_metadata.dt <- fread(io$metacell_metadata)
# sample_metadata.dt <- fread(io$metadata)

###############################
## Load SingleCellExperiment ##
###############################

sce <- readRDS(io$metacell_sce)

###########################
## Load precomputed UMAP ##
###########################

# umap.dt <- fread(io$umap, select=c(3,1,2)) %>% setnames(c("cell","V1","V2"))

umap.dt <- fread(io$rna.atlas.metadata) %>%
  .[stripped==F & doublet==F] %>%
  .[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

#########################################################
## Plot dimensionality reduction: one sample at a time ##
#########################################################

to.plot <- umap.dt %>% copy %>%
  .[,index:=match(cell, metacell_metadata.dt$closest.cell)] %>% 
  .[,mapped:=as.factor(!is.na(index))] %>% 
  .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas","Metacell"))] %>%
  setorder(mapped) 

p <- plot.dimred(to.plot, query.label = "Metacell", atlas.label = "Atlas")

pdf(file.path(io$outdir,"umap_metacell.pdf"), width=8, height=6.5)
print(p)
dev.off()
