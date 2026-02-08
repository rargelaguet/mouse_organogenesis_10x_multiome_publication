library(scater)
library(scran)
library(purrr)
library(ggplot2)
library(data.table)

sce <- readRDS("/home/irosshandler/partial-reprogramming/human/lung/HumanLungCellAtlas_Nature/sce_rawcounts.rds")


head(colData(sce))
colnames(colData(sce))

table(sce$ann_level_1)

sce
counts(sce)[1:5,1:5]

unique(sce$cell_type)

# Subset fibroblasts
sce <- sce[,sce$ann_finest_level%in% c("Adventitial fibroblasts", "Alveolar fibroblasts","Peribronchial fibroblasts","Myofibroblasts","Fibromyocytes")]

# Filter donors with less than 100 cells
tmp <- table(sce$donor_id)
donors.to.use <- names(which(tmp>=100))
sce <- sce[,sce$donor_id%in%donors.to.use]

table(sce$age_or_mean_of_age_range)
table(sce$cell_type)
table(sce$dataset)

sort(table(sce$ann_finest_level,sce$cell_type)["Adventitial fibroblasts",])
sort(table(sce$ann_finest_level,sce$cell_type)["Alveolar fibroblasts",])
sort(table(sce$ann_finest_level,sce$cell_type)["Peribronchial fibroblasts",])
sort(table(sce$ann_finest_level,sce$cell_type)["Peribronchial fibroblasts",])
sort(table(sce$ann_finest_level,sce$cell_type)["Myofibroblasts",])
sort(table(sce$ann_finest_level,sce$cell_type)["Fibromyocytes",])

saveRDS(sce, "/share/rargelaguet/data/Sikkema2023_lung_cell_atlas/processed/SingleCellExperiment_fibroblasts.rds")


#############
## Explore ##
#############

unique(sce$donor_id[sce$age_or_mean_of_age_range>=60])
dput(unique(sce$donor_id))

# Subset
sce.tmp <- sce[,sce$donor_id%in%"homosapiens_None_2023_None_sikkemalisa_002_d10_1101_2022_03_10_483747A44"]
sce.tmp <- sce[,sce$donor_id%in%"homosapiens_None_2023_None_sikkemalisa_002_d10_1101_2022_03_10_483747VUHD94"]

# Normalisation
sce.tmp <- scuttle::logNormCounts(sce.tmp)

# Feature selection
decomp <- modelGeneVar(sce.tmp)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=2500) %>% rownames
sce_filt.tmp <- sce.tmp[hvgs,]

# PCA
sce_filt.tmp <- runPCA(sce_filt.tmp, ncomponents=2, ntop=2500)
reducedDims(sce.tmp) <- reducedDims(sce_filt.tmp)

# Plot
plotReducedDim(sce.tmp, dimred = "PCA", colour_by="cell_type", point_size=2, point_alpha=1)
plotReducedDim(sce.tmp, dimred = "PCA", colour_by="ENSG00000111640", point_size=2, point_alpha=1) # GAPDH
plotReducedDim(sce.tmp, dimred = "PCA", colour_by="ENSG00000026025", point_size=2, point_alpha=1) # VIM
plotReducedDim(sce.tmp, dimred = "PCA", colour_by="ENSG00000245532", point_size=2, point_alpha=1) # NEAT1 
plotReducedDim(sce.tmp, dimred = "PCA", colour_by="ENSG00000135404", point_size=2, point_alpha=1) # CD63
plotReducedDim(sce.tmp, dimred = "PCA", colour_by="ENSG00000148773", point_size=2, point_alpha=1) # MKI67
plotReducedDim(sce.tmp, dimred = "PCA", colour_by="ENSG00000108821", point_size=2, point_alpha=1) # COL1A1
plotReducedDim(sce.tmp, dimred = "PCA", colour_by="ENSG00000160789", point_size=2, point_alpha=1) # LMNA

plot()
