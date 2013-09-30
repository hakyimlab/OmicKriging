\name{load_gene_data}
\alias{load_gene_data}
\title{Load genetic data from PLINK binary files.}
\usage{
  load_gene_data(bedFile, bimFile, famFile, gdsFile)
}
\arguments{
  \item{bedFile}{File path to PLINK .bed file.}

  \item{bimFile}{File path to PLINK .bim file.}

  \item{famFile}{File path to PLINK .fam file.}

  \item{gdsFile}{File path to store GDS for other
  analyses.}
}
\value{
  None. The GDS file is stored on disk.
}
\description{
  This function takes file paths for PLINK binary files and
  loads them into a Genomic Data Structure (GDS) stored on
  disk. This is a wrapper function aroud utilities provided
  by the gdsfmt and SNPRelate packages. The GDS is used in
  downstream analyses such as: computing the genetic
  related matrix (\code{\link{make_grm}}) and computing
  principal components (\code{\link{make_PCs}})
}
\references{
  library(gdsfmt), library(SNPRelate)
}
\keyword{input}

