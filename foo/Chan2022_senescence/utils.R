# load_SummarizedExperiment <- function(file, normalise = FALSE, features = NULL, cells = NULL, remove_non_expressed_genes = FALSE) {
#   sce <- readRDS(file)
#   if (!is.null(cells)) sce <- sce[,cells]
#   if (!is.null(features)) sce <- sce[features,]
#   if (normalise) sce <- scuttle::logNormCounts(sce)
#   if (remove_non_expressed_genes) sce <- sce[which(Matrix::rowMeans(counts(sce))>1e-4),]
#   return(sce)
# }

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}