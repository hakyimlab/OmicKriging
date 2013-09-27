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
make_PCs <- function(gdsFile, n.core, n.top = 0) {
  gds <- openfn.gds(gdsFile)
  pca <- snpgdsPCA(gds, num.thread = n.core)
  return( pca$eigenvect )
}
