#' Compute genetic correlation matrix from PLINK binary files.
#'
#' This is a convenience function which produces a centered genetic correlation
#' matrix from SNPs in PLINK binary files. The resulting matrix can be used 
#' with the \code{\link{okriging}} function. It uses the Genomic Data Structure
#' for fast I/O.
#'
#' @param bedFile PLINK bed file.
#' @param bimFile PLINK bim file.
#' @param famFile PLINK fam file.
#' @param gdsFile File to store the Genomic Data Structure on disk. Default is
#'   a temporary file created by tempfile().
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
make_grm <- function(bedFile, bimFile, famFile, gdsFile = tempfile(), grmDataFile = NULL) {
  require(gdsfmt)
  require(SNPRelate)
  source('R/rcppcormat.r')

  ## convert binary files to GDS file
  snpgdsBED2GDS(bedFile, famFile, bimFile, gdsFile)
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
