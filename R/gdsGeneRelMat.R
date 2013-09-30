#' Compute genetic correlation matrix from PLINK binary files using the Genomic Data Structure (GDS).
#'
#' This is a convenience function which produces a centered genetic correlation
#' matrix from SNPs loaded into a Genomic Data Structure (GDS) file. This function differs
#' from \code{\link{make_grm}} in that it leverages the gdsfmt library for the
#' parallel calculation of the covariance matrix. The resulting matrix can be used 
#' with the \code{\link{okriging}} function. The GRM can be saved to disk as a
#' R object for fast loading downstream.
#'
#' @param gdsFile File to store the Genomic Data Structure on disk for use elsewhere.
#' @param grmDataFile File to store the resulting GRM on disk as an R object.
#' @param snpList A vector of SNP IDs to subset the GRM on.
#' @param sampleList A vector of sample IDs to subset the GRM on.
#' @param n.core Number of cores to distribute computation across.
#'
#' @return A genetic correlation matrix with colnames and rownames set to sample IDs.
#'   Each entry in the matrix is of type 'double'.
#'
#' @import gdsfmt
#' @import SNPRelate
#'
#' @keywords input, GRM
#' @export
make_grm_gds <- function(gdsFile, grmDataFile = NULL, sampleList = NULL, snpList = NULL, n.core = 1) {
  require(gdsfmt)
  require(SNPRelate)

  genofile <- openfn.gds(gdsFile)
  ## compute the genetic covariance matrix
  gcov <- snpgdsPCA(genofile, sample.id = sampleList, snp.id = snpList, num.thread = n.core, genmat.only = TRUE, verbose = FALSE)
  ## make correlation matrix from covariance matrix
  gcor <- cov2cor(gcov$genmat)

  ## save to file if provided
  if(!is.null(grmDataFile)) {
    save(gcor, file = grmDataFile)
  }

  ## pull sample IDs unless a sample list is specified
  if(!is.null(sampleList)) {
    sample.ids <- sampleList
  } else {
    sample.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
  }
  ## add sample IDs to GRM object
  rownames(gcor) <-  sample.ids
  colnames(gcor) <-  sample.ids

  return(gcor)
}

