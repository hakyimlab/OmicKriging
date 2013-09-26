## function to compute the genetic relatedness matrix from plink binary files

## *******************************
rcppcormat <- function(snpmat){
    ## require
    require(Rcpp)
    require(RcppEigen)
    require(inline)
    ## rcpp
    crossprodCpp <- '
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::Lower;

    // set constants

    const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
    const int n(A.cols());

    // create n x n matrix

    MatrixXd AtA(MatrixXd(n, n)

    // zero out the matrix

    .setZero()

    // treat what follows as symmetric matrix, use lower

    .selfadjointView<Lower>()
    // sum of B + AA

    .rankUpdate(A.adjoint()));

    // return

    return wrap(AtA);
    '

    ## compile the cross product function
    cpcpp <- cxxfunction(signature(AA="matrix"), crossprodCpp,
    plugin="RcppEigen", verbose=FALSE)

    cormat <- cpcpp(snpmat)

    ## post work
    cormatdiag <- diag(cormat)
    cormat <- sweep(cormat, 1, cormatdiag, "/")
    cormat <- sweep(cormat, 2, cormatdiag, "/")
    return(cormat)

}

## compute genetic correlation matrix from plink binary files
make_grm <- function(bedFile, bimFile, famFile, gdsFile=tempfile()) {

  require(gdsfmt)
  require(SNPRelate)

  ## convert binary files to GDS file
  snpgdsBED2GDS(bedFile, famFile, bimFile, gdsFile)
  genofile <- openfn.gds(gdsFile)

  ## pull sample IDs (SNP-major mode)
  sample.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))

  ## pull genotype matrix, center, and compute the resulting correlation matrix
  X <- read.gdsn(index.gdsn(genofile, "genotype"))
  Xbar <- sweep(X, 2, colMeans(X), "-")
  grm <- rcppcormat(ti(Xbar))

  ## annotate columns and rows with sample IDs
  colnames(grm) <- sample.ids
  rownames(grm) <- sample.ids

  return(grm)
}
