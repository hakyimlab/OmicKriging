#' Compute genetic correlation matrix from PLINK binary files.
#'
#' This is a convenience function which produces a centered genetic correlation
#' matrix from SNPs loaded into a Genomic Data Structure (GDS) file. The resulting matrix can be used 
#' with the \code{\link{okriging}} function. The GRM can be saved to disk as a
#' R object for fast loading downstream. The calculation first subtracts allele
#' dosage (i.e. column) means from each field, and sets missing values to column
#' mean (i.e. 0).
#'
#' @param gdsFile File to store the Genomic Data Structure on disk for use elsewhere.
#' @param grmDataFile File to store the resulting GRM on disk as an R object.
#' @param snpList A vector of SNP IDs to subset the GRM on.
#' @param sampleList A vector of sample IDs to subset the GRM on.
#'
#' @return A genetic correlation matrix with colnames and rownames set to sample IDs.
#'   Each entry in the matrix is of type 'double'.
#'
#' @include R/rcppcormat.r
#'
#' @import gdsfmt
#' @import SNPRelate
#'
#' @keywords input, GRM
#' @examples
#' grm <- make_grm(bedFile = "data/T1DCC.subset.bed",
#'          bimFile = "data/T1DCC.subset.bim",
#'          famFile = "data/T1DCC.subset.fam",
#'          gdsFile = "~/tmp/T1DCC.subset.gds")
#' @export
make_grm <- function(gdsFile = NULL, grmDataFile = NULL, snpList = NULL, sampleList = NULL) {
  require(gdsfmt)
  require(SNPRelate)
  source('R/rcppcormat.r')

  genofile <- openfn.gds(gdsFile)
  ## pull an integer dosage matrix from the GDS. Rows are samples, columns are SNPs, and missing values are int 3.
  X <- snpgdsGetGeno(gdsobj = genofile, sample.id = sampleList, snp.id = snpList, verbose = FALSE)
  ## set missing values (int 3) to properly missing
  X[X == 3] <- NA
  Xbar <- sweep(X, 2, colMeans(X, na.rm = TRUE), "-")
  X[X == 3] <- 0.0
  grm <- rcppcormat(t(Xbar))
  
  ## pull sample IDs unless a sample list is specified
  if(is.null(sampleList)) {
    sample.ids <- sampleList
  } else {
    sample.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
  }

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
