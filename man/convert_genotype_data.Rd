\name{convert_genotype_data}
\alias{convert_genotype_data}
\title{Load genetic data from PLINK binary files.}
\usage{
  convert_genotype_data(bFile, gdsFile)
}
\arguments{
  \item{bFile}{File path to PLINK .bed .bim .fam files.}

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
  related matrix (make_GRM) and computing principal
  components (make_PCA).
}
\references{
  library(gdsfmt), library(SNPRelate)
}
\keyword{input}

