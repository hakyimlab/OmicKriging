\name{read_GRMBin}
\alias{read_GRMBin}
\title{Read the GRM binary file.}
\usage{
  read_GRMBin(prefix, size = 4)
}
\arguments{
  \item{prefix}{The file path prefix to GRM binary files
  (e.g., test.grm.bin, test.grm.N.bin, test.grm.id.)}

  \item{size}{The length (in bytes) of each value in the
  raw GRM vector. Default is 4, and matches GRM writen by
  GCTA 1.11.}
}
\value{
  GRM of dim (N.samples x N.samples) with rownames and
  colnames as sample ID.
}
\description{
  Function provided by GCTA maintainers (modified slightly)
  for accessing their recently introduced binary GRM
  format. The GRM is stored as a vector of numerics which
  correspond to the lower triangular elements including the
  diagonal. We simply read these, pull the diagonal
  elements, and inflate them into a full symmetric matrix.
  We add sample IDs to colnames and rownames for
  compatibility with other Kriging functions.
}
\references{
  http://www.complextraitgenomics.com/software/gcta/estimate_grm.html
}

