champ.se <- readRDS("/Users/rargelaguet/data/Chan2022_senescence/processed/atac/SummarizedExperiment_ATAC_Jack.rds") 
rowData(champ.se)
colData(champ.se)

reducedDim(champ.se)

assayNames(champ.se)
hist(as.matrix(assay(champ.se,"rawCounts")))
hist(as.matrix(assay(champ.se,"normCounts")))



sample_metadata.dt <- colData(champ.se) %>% as.data.table(keep.rownames=T) %>%
  setnames("rn","sample")
fwrite(sample_metadata.dt, "/Users/rargelaguet/data/Chan2022_senescence/processed/atac/sample_metadata.txt.gz", quote=F, sep="\t", na="NA")


feature_metadata.dt <- rowData(champ.se) %>% as.data.table(keep.rownames=F) %>%
  .[,c("HOMER_NearestGeneDesc","HOMER_DetailedAnnotations"):=NULL]
fwrite(feature_metadata.dt, "/Users/rargelaguet/data/Chan2022_senescence/processed/atac/feature_metadata.txt.gz", quote=F, sep="\t", na="NA")
