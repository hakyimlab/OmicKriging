#' Compute genetic correlation matrix from PLINK binary files.
#'
#' This is a convenience function which produces a centered genetic correlation
#' matrix from SNPs loaded into a Genomic Data Structure (GDS) file. The resulting matrix can be used 
#' with the \code{\link{okriging}} function. The GRM can be saved to disk as a
#' R object for fast loading downstream.
#'
#' @param gdsFile File to store the Genomic Data Structure on disk for use elsewhere.
#' @parma grmDataFile File to store the resulting GRM on disk as an R object.
#'
#' @return A genetic correlation matrix with colnames and rownames set to sample IDs.
#'
#' @include R/rcppcormat.r
#'
#' @keywords input
#' @examples
#' grm <- make_grm(bedFile = "data/T1DCC.subset.bed",
#'          bimFile = "data/T1DCC.subset.bim",
#'          famFile = "data/T1DCC.subset.fam",
#'          gdsFile = "~/tmp/T1DCC.subset.gds")
make_grm <- function(gdsFile = NULL, grmDataFile = NULL) {
  require(gdsfmt)
  require(SNPRelate)
  source('R/rcppcormat.r')

  genofile <- openfn.gds(gdsFile)

  ## pull sample IDs (SNP-major mode)
  sample.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))

  ## pull genotype matrix, center, and compute the resulting correlation matrix
  X <- read.gdsn(index.gdsn(genofile, "genotype"))
  Xbar <- sweep(X, 2, colMeans(X), "-")
  grm <- rcppcormat(t(Xbar))

  ## annotate columns and rows with sample IDs
  colnames(grm) <- sample.ids
  rownames(grm) <- sample.ids

  ## write out the GRM if a file is specified
  if( !is.null(grmDataFile) ) {
    save(grm, file = grmDataFile)
  }

  return(grm)
}

#' Load genetic relatedness matrix from file.
#'
#' Loads a genetic relatedness matrix (GRM) from an .Rdata file produced by save()
#' in the \code{\link{make_grm}} function. This can be used to quickly load a 
#' GRM instead of recomputing it from PLINK binary files
#'
#' @param grmDataFile Location of the .Rdata file.
#'
#' @return A genetic correlation matrix with colnames and rownames set to sample IDs.
#'
#' @keywords input
#' 
#' @export
load_grm <- function(grmDataFile) {
  return(load(grmDataFile))
}
