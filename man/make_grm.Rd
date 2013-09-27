\name{make_grm}
\alias{make_grm}
\title{Compute genetic correlation matrix from PLINK binary files.}
\usage{
  make_grm(bedFile, bimFile, famFile, gdsFile = tempfile())
}
\arguments{
  \item{bedFile}{PLINK bed file.}

  \item{bimFile}{PLINK bim file.}

  \item{famFile}{PLINK fam file.}

  \item{gdsFile}{File to store the Genomic Data Structure
  on disk.}
}
\value{
  A genetic correlation matrix with colnames and rownames
  set to sample IDs.
}
\description{
  This is a convenience function which produces a centered
  genetic correlation matrix from SNPs in PLINK binary
  files. The resulting matrix can be used with the
  \code{\link{okriging}} function. It uses the Genomic Data
  Structure for fast I/O.
}
\examples{
grm <- make_grm(bedFile = "data/T1DCC.subset.bed",
         bimFile = "data/T1DCC.subset.bim",
         famFile = "data/T1DCC.subset.fam",
         gdsFile = "~/tmp/T1DCC.subset.gds")
}
\keyword{input}

