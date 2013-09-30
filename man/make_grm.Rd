\name{make_grm}
\alias{make_grm}
\title{Compute genetic correlation matrix from PLINK binary files.}
\usage{
  make_grm(gdsFile = NULL, grmDataFile = NULL,
    snpList = NULL, sampleList = NULL)
}
\arguments{
  \item{gdsFile}{File to store the Genomic Data Structure
  on disk for use elsewhere.}

  \item{grmDataFile}{File to store the resulting GRM on
  disk as an R object.}

  \item{snpList}{A vector of SNP IDs to subset the GRM on.}

  \item{sampleList}{A vector of sample IDs to subset the
  GRM on.}
}
\value{
  A genetic correlation matrix with colnames and rownames
  set to sample IDs.  Each entry in the matrix is of type
  'double'.
}
\description{
  This is a convenience function which produces a centered
  genetic correlation matrix from SNPs loaded into a
  Genomic Data Structure (GDS) file. The resulting matrix
  can be used with the \code{\link{okriging}} function. The
  GRM can be saved to disk as a R object for fast loading
  downstream.
}
\examples{
grm <- make_grm(bedFile = "data/T1DCC.subset.bed",
         bimFile = "data/T1DCC.subset.bim",
         famFile = "data/T1DCC.subset.fam",
         gdsFile = "~/tmp/T1DCC.subset.gds")
}
\keyword{GRM}
\keyword{input,}

