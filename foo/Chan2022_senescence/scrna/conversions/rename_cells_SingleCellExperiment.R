library(SingleCellExperiment)

cell_metadata.dt <- fread("/Users/argelagr/data/Xing2020_bj_oskm/cell_metadata.txt.gz")

sce <- readRDS("/Users/argelagr/data/Xing2020_bj_oskm/rna/processed/SingleCellExperiment.rds")
tmp <- cell_metadata.dt[cell_old%in%colnames(sce)] %>% setkey(cell_old) %>% .[colnames(sce)]
stopifnot(colnames(sce)==tmp$cell_old)
colnames(sce) <- tmp$cell

saveRDS(sce, "/Users/argelagr/data/Xing2020_bj_oskm/rna/processed/SingleCellExperiment.rds")
