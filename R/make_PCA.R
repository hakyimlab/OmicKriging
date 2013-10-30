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
#' @import SNPRelate
#'
#' @references library(SNPRelate)
#' 
#' @export
make_PCs_gds <- function(gdsFile, n.core, n.top = 0) {

  gds <- openfn.gds(gdsFile)
  pca <- snpgdsPCA(gds, num.thread = n.core)
  rownames(pca$eigenvect) <- read.gdsn(index.gdsn(gds, "sample.id"))
  return( pca$eigenvect )
}

#' Run Principal Component Analysis (PCA) using base R svd() function.
#'
#' A simple wrapper around the base R svd() function which returns the top N
#' eigenvectors of a matrix. Use this function to generate covariates for use
#' with the \code{\link{okriging}} or \code{\link{krigr_cross_validation}}
#' functions.
#'
#' @param X A correlation matrix.
#' @param n.top Number of top principal compenents to return 
#'
#' @return A matrix of Principal Components of dimension (# of samples) x
#'   (n.top). As expected, eigenvectors are ordered by eigenvalue. Rownames
#'   are given as sample IDs.
#'
#' @keywords covariate, PCA, GRM
#' @export
make_PCs_svd <- function(X, n.top = 2) {
  res <- La.svd(X, nu = n.top)
  rownames(res) <- rownames(X)
  return(res["u"])
}

#' Run Principal Component Analysis (PCA) using the irlba package.
#'
#' A simple wrapper around the irlba() function which computes a partial SVD
#' efficiently. This function's run time depends on the number of eigenvectors
#' requested but scales well. Use this function to generate covariates for use
#' with the \code{\link{okriging}} or \code{\link{krigr_cross_validation}}
#' functions.
#'
#' @param X A correlation matrix.
#' @param n.top Number of top principal compenents to return 
#'
#' @return A matrix of Principal Components of dimension (# of samples) x
#'   (n.top). As expected, eigenvectors are ordered by eigenvalue. Rownames
#'   are given as sample IDs.
#'
#' @keywords covariate, PCA, GRM
#'
#' @references library(irlba)
#'
#' @import irlba
#' @export
make_PCs_irlba <- function(X, n.top = 2) {
  
  res <- irlba(X, nu = n.top)
  rownames(res$u) <- rownames(X)
  return(res$u)
}

