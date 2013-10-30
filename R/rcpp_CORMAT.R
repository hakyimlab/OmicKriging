#' Compute a correlation matrix.
#'
#' This function computes a correlation matrix using a cross product. It is
#' implemented in C++ for performance. This function is compiled and called
#' from R using the Rcpp package.
#'
#' @param snpmat A centered SNP matrix containing fields of type 'double'.
#'
#' @return A correlation matrix.
#'
#' @keywords c++, performance, correlation matrix
rcppcormat <- function(snpmat){
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
