## function to compute the genetic relatedness matrix using the SNPRelate package

make_grm_gds <- function(gdsFile, grmDataFile = NULL, sampleList = NULL, snpList = NULL, n.core = 1) {
  require(gdsfmt)
  require(SNPRelate)

  genofile <- openfn.gds(gdsFile)
  ## compute the genetic covariance matrix
  gcov <- snpgdsPCA(genofile, sample.id = sampleList, snp.id = snpList, num.thread = n.core, genmat.only = TRUE, verbose = FALSE)
  ## make correlation matrix from covariance matrix
  gcor <- cov2cor(gcov$genmat)

  ## save to file if provided
  if(!is.null(grmDataFile)) {
    save(gcor, file = grmDataFile)
  }

  ## pull sample IDs unless a sample list is specified
  ed.gdsn(index.gdsn(genofile, "sample.id"))f(!is.null(sampleList)) {
    sample.ids <- sampleList
  } else {
    sample.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
  }
  ## add sample IDs to GRM object
  rownames(gcor) <-  sample.ids
  colnames(gcor) <-  sample.ids

  return(gcor)
}

