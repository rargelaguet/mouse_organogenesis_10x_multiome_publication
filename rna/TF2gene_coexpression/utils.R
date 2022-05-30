library(Matrix)
library(qlcMatrix)
set.seed(42)

#' Replace non-zero entries in a sparse entries with non-zero ranks
#'
#' This method creates a rank matrix for a sparse matrix X using the following approach:
#' 1. Use non-zero enries in a column to calculate the ranks
#' 2. Add (z-1)/2 to the ranks (only non-zero entries are changed). z is the number of zeros
#' in the column
#' Since all the entries are shifted by the same constant (the zeros
#' are already shifted), the covariance matrix of this shifted matrix is
#' the same as the rank matrix of the entire matrix (where the zeros would
#' all also have a rank = (z+1)/2) where z is the number of zeros
#'
#' This rank matrix can then be used to calculate pearson correlation

SparsifiedRanks <- function(X) {
  X <- as(object = X, Class = "dgCMatrix")
  j <- summary(object = X)$j
  n_zeros_per_col <- nrow(X) - diff(X@p)

  for (column in unique(x = j)) {
    non_zero_element_index <- which(j == column)
    elements_along_row <- X@x[non_zero_element_index]
    ranks <- rank(elements_along_row)
    ranks <- ranks + (n_zeros_per_col[column] - 1) / 2
    X@x[non_zero_element_index] <- ranks
  }
  return(X)
}

SparseSpearmanCor <- function(X, Y = NULL, cov = FALSE) {

  # Get sparsified ranks
  rankX <- SparsifiedRanks(X)
  if (is.null(Y)){
    # Calculate pearson correlation on rank matrices
    return (corSparse(X=rankX, cov=cov))
    }
  rankY <- SparsifiedRanks(Y)
  return(corSparse( X = rankX, Y = rankY, cov = cov))
}


SparsifiedRanks2 <- function(X) {
  if (class(X)[1] != "dgCMatrix") {
    X <- as(object = X, Class = "dgCMatrix")
  }
  non_zeros_per_col <- diff(x = X@p)
  n_zeros_per_col <- nrow(x = X) - non_zeros_per_col
  offsets <- (n_zeros_per_col - 1) / 2
  x <- X@x
  ## split entries to columns
  col_lst <- split(x = x, f = rep.int(1:ncol(X), non_zeros_per_col))
  ## calculate sparsified ranks and do shifting
  sparsified_ranks <- unlist(x = lapply(X = seq_along(col_lst), FUN = function(i) rank(x = col_lst[[i]]) + offsets[i]))
  ## Create template rank matrix
  X.ranks <- X
  X.ranks@x <- sparsified_ranks
  return(X.ranks)
}


SparseSpearmanCor2 <- function(X, Y = NULL, cov = FALSE) {

  # Get sparsified ranks
  rankX <- SparsifiedRanks2(X)
  if (is.null(Y)){
    # Calculate pearson correlation on rank matrices
    return (corSparse(X=rankX, cov=cov))
    }
  rankY <- SparsifiedRanks2(Y)
  return(corSparse( X = rankX, Y = rankY, cov = cov))
}