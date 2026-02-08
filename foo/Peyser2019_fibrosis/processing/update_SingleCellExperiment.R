library(scater)
library(scran)
library(purrr)
library(ggplot2)
library(data.table)


sce <- readRDS("/Users/rargelaguet/data/Sikkema2023_lung_cell_atlas/processed/SingleCellExperiment_fibroblasts.rds")

head(colData(sce))
colnames(colData(sce))


###################
## Cell metadata ##
###################

#######################
## Update gene names ##
#######################

gene_metadata.dt <- fread("/Users/rargelaguet/data/ensembl/human/v99/bioMart/all_genes/Hsapiens_genes_BioMart.99.txt.gz") %>% .[!duplicated(symbol)]

#######################
## Update gene names ##
#######################

ens_id.to.use <- intersect(gene_metadata.dt$ens_id,rownames(sce))
gene_metadata.dt <- gene_metadata.dt[ens_id%in%ens_id.to.use,c("ens_id","symbol")] %>% setkey(ens_id) %>% .[ens_id.to.use]
sce <- sce[ens_id.to.use,]
stopifnot(rownames(sce)==gene_metadata.dt$ens_id)
rownames(sce) <- gene_metadata.dt$symbol
tmp <- gene_metadata.dt %>% DataFrame; rownames(tmp) <- tmp$symbol
rowData(sce) <- tmp

##############################
## Dimensionality reduction ##
##############################

reducedDimNames(sce)

plotReducedDim(sce, dimred = "UMAP", colour_by="GAPDH", point_size=2, point_alpha=1, by_exprs_values="counts")

##########
## Save ##
##########

cell_metadata.dt <- colData(sce) %>% as.data.table(keep.rownames=F)# %>% setnames("rn","cell")

fwrite(cell_metadata.dt, "/Users/rargelaguet/data/Sikkema2023_lung_cell_atlas/processed/cell_metadata.txt.gz", sep="\t", quote=T, na="NA")
saveRDS(sce, "/Users/rargelaguet/data/Sikkema2023_lung_cell_atlas/processed/SingleCellExperiment_fibroblasts.rds")