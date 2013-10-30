\name{make_GRM}
\alias{make_GRM}
\title{Compute genetic correlation matrix from PLINK binary files.}
\usage{
  make_GRM(gdsFile = NULL, grmFilePrefix = NULL,
    snpList = NULL, sampleList = NULL)
}
\arguments{
  \item{gdsFile}{File holding the GDS from which to pull
  the raw genotype matrix.}

  \item{grmFilePrefix}{File to store the resulting GRM on
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
  downstream. The genotype matrix is z-score normalized
  (i.e. column means are centered and column variance is
  divided out to unit variance) prior to calculating the
  correlation matrix.
}
\keyword{GRM}
\keyword{input,}

