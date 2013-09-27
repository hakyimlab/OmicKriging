## function to compute the genetic relatedness matrix from plink binary files


## compute genetic correlation matrix from plink binary files
make_grm <- function(bedFile, bimFile, famFile, gdsFile=tempfile()) {

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

  return(grm)
}
