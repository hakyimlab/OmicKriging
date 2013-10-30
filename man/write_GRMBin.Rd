\name{write_GRMBin}
\alias{write_GRMBin}
\title{Write GRM binary files.}
\usage{
  write_GRMBin(X, n.snps = 0, prefix, size = 4)
}
\arguments{
  \item{X}{Correlation matrix with rownames and colnames as
  sample IDs.}

  \item{prefix}{Base file path and names for the three
  output files.}

  \item{n.snps}{Number of SNPs used in correlation matrix
  calculation. Default is 0.0.}

  \item{size}{Number of bytes to write for each value.
  Default is 4}
}
\value{
  None. Though side effects are writing three files as
  described above.
}
\description{
  Function to write a binary GRM format recently introduced
  by GCTA. It takes a correlation matrix as used by other
  Kriging functions, and writes three files: binary file
  for storing the diagonal + lower triangular elements, a
  text file for sample IDs, and a binary file storing the
  number of SNPs used in the correlation matrix
  calculation.
}
\references{
  http://www.complextraitgenomics.com/software/gcta/estimate_grm.html
}

