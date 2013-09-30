#' Run Principal Component Analysis (PCA) using the Genomic Data Structure (GDS).
#'
#' An efficient method for computing Principal Components using the Genomic Data
#' Structure (GDS). This is a convenience wrapper for functions from the 
#' SNPRelate package.
#'
#' @param gdsFile A Genomic Data Structure file describing your study.
#' @param n.top Number of top principal components to return. Defaults to
#'   returning all components (i.e. # of samples).
#' @param n.core Distrubute computation across N cores.
#'
#' @return A matrix of Principal Components of dimension (# of samples) x
#'   (n.top).
#'
#' @keywords covariate, PCA, GRM
#'
#' @reference library(SNPRelate)
#' 
#' @export
make_PCs_gds <- function(gdsFile, n.core, n.top = 0) {
  gds <- openfn.gds(gdsFile)
  pca <- snpgdsPCA(gds, num.thread = n.core)
  return( pca$eigenvect )
}

## SVD based PCA. Assuming that the matrix X is a valid correlation matrix with rownames as sample names.
make_PCs_svd <- function(X, n.top = 2) {
  res <- La.svd(X, nu = n.top)
  rownames(res) <- rownames(X)
  return(res["u"])
}

## IRLBA based SVD -- supposedly the state of the art.
make_PCs_irlba <- function(X, n.top = 2) {
  require(irlba)
  
  res <- irlba(X, nu = n.top)
  rownames(res) <- rownames(X)
  return(res["u"])
}

