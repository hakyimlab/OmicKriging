#' Compute genetic correlation matrix from PLINK binary files.
#'
#' This is a convenience function which produces a centered genetic correlation
#' matrix from SNPs loaded into a Genomic Data Structure (GDS) file. The resulting matrix can be used 
#' with the okriging function. The GRM can be saved to disk as a
#' R object for fast loading downstream. The genotype  matrix is z-score normalized (i.e.
#' column means are centered and column variance is divided out to unit variance)
#' prior to calculating the correlation matrix.
#'
#' @param gdsFile File holding the GDS from which to pull the raw genotype matrix.
#' @param grmFilePrefix File to store the resulting GRM on disk as an R object.
#' @param snpList A vector of SNP IDs to subset the GRM on.
#' @param sampleList A vector of sample IDs to subset the GRM on.
#'
#' @return A genetic correlation matrix with colnames and rownames set to sample IDs.
#'   Each entry in the matrix is of type 'double'.
#'
#' @include R/grm_IO.R
#' @include R/rcpp_CORMAT.R
#'
#' @import gdsfmt
#' @import SNPRelate
#'
#' @keywords input, GRM
#' @export
make_GRM <- function(gdsFile = NULL, grmFilePrefix = NULL, snpList = NULL, sampleList = NULL) {
  source('R/grm_IO.R')
  source('R/rcpp_CORMAT.R')

  genofile <- openfn.gds(gdsFile)
  ## pull an integer dosage matrix from the GDS. Rows are samples, columns are SNPs, and missing values are int 3.
  X <- snpgdsGetGeno(gdsobj = genofile, sample.id = sampleList, snp.id = snpList, verbose = FALSE)
  ## set missing values (int 3) to properly missing
  X[X == 3] <- NA
  ## z-normalize matrix (sweep out column means, and divide out column matrices)
  X <- scale(X, center = TRUE, scale = TRUE)
  ## set missing values to new column mean, i.e. 0.0
  X[is.na(X)] <- 0.0
  
  gc()
  
  grm <- rcppcormat(t(X))
  
  ## pull sample IDs unless a sample list is specified
  if(!is.null(sampleList)) {
    sample.ids <- sampleList
  } else {
    sample.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
  }

  ## annotate columns and rows with sample IDs
  colnames(grm) <- sample.ids
  rownames(grm) <- sample.ids

  ## write out the GRM if a file is specified
  if(!is.null(grmFilePrefix)) {
    write_GRMBin(X = grm, prefix = grmFilePrefix)
  }

  return(grm)
}
