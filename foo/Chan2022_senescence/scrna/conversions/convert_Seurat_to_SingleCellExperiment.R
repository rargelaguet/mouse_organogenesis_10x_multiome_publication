library(Seurat)
library(SingleCellExperiment)

seurat <- readRDS("/Users/rargelaguet/data/Xing2020_bj_oskm/processed/rna/seurat.rds")

cell_metadata.dt <- fread("/Users/rargelaguet/data/Xing2020_bj_oskm/cell_metadata.txt.gz")
cell_metadata.dt[!cell_old%in%colnames(seurat)]

sce <- as.SingleCellExperiment(seurat)
logcounts(sce) <- NULL
assays(sce)[["scaledata"]] <- NULL

head(colData(sce))
head(rowData(sce))

head(colnames(sce))
head(rownames(sce))

# cell_metadata.dt <- fread("/Users/argelagr/data/Xing2020_bj_oskm/raw_data/cell_metadata.txt.gz") %>%
cell_metadata_filt.dt <- cell_metadata.dt %>%
  .[assay=="snRNA"] %>%
  .[,sample_id:=strsplit(cell,"_") %>% map_chr(1)] %>%
  .[grepl("-FM|-PR",sample_id)]

foo <- cell_metadata.dt$sample_id %>% unique
bar <- strsplit(colnames(sce),"_") %>% map_chr(1) %>% unique

colnames(sce) <- gsub("D0_","D0-FM_",colnames(sce))
colnames(sce) <- gsub("D4_","D4-FM_",colnames(sce))
colnames(sce) <- gsub("D8_","D8-FM_",colnames(sce))
colnames(sce) <- gsub("D12_","D12-PR_",colnames(sce))
colnames(sce) <- gsub("D16_","D16-PR_",colnames(sce))
colnames(sce) <- gsub("D20_","D20-PR_",colnames(sce))
colnames(sce) <- gsub("D24_","D24-PR_",colnames(sce))
colnames(sce) <- gsub("P20_","P20-PR_",colnames(sce))

cells <- intersect(colnames(sce),cell_metadata.dt$cell)
# table(cell_metadata.dt[!cell%in%cells]$sample_id)
cell_metadata_to_sce.dt <- cell_metadata.dt[cell%in%cells] %>% setkey(cell) %>% .[cells]
sce <- sce[,cells]
colData(sce) <- cell_metadata_to_sce.dt %>% tibble::column_to_rownames("cell") %>% DataFrame


saveRDS(sce, "/Users/argelagr/data/Xing2020_bj_oskm/processed/rna/SingleCellExperiment.rds")
