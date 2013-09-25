computeGX <- function(genefile, outheader, idfile=NULL){
    ## functions
    require(Rcpp)
    require(RcppEigen)
    require(inline)
    
    ## cross product
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

    // sum of B + AAt

    .rankUpdate(A.adjoint()));

    // return

    return wrap(AtA);
    '

    ## compile the cross product function
    cpcpp <- cxxfunction(signature(AA="matrix"), crossprodCpp,
    plugin="RcppEigen", verbose=FALSE)
    
    ##
    genedata <- read.delim(genefile, sep="", as.is=T, header=T)
    if(!is.null(idfile))
    {
    iddata <- read.table(idfile,header=T,as.is=T)
    genedata <- merge(iddata,genedata,by.x=c("FID","IID"),by.y=c("FID","IID"))
    }
    cordata.id <- genedata[c("FID","IID")]
    genemat <- as.matrix(genedata[,!(names(genedata) %in% c("FID","IID"))])

    ## impute missing data
    if(sum(is.na(genemat))>0) print("Impute first") ## TODO: impute with the script
    cormat <- cpcpp(t(genemat))
    cormatdiag <- diag(cormat)
    cormat <- sweep(cormat, 1, cormatdiag, "/")
    cormat <- sweep(cormat, 2, cormatdiag, "/")

    ## format matrix into linear dataframe (GCTA format)
    ng <- ncol(genemat)
    cordata <- cbind(c(col(cormat)),c(row(cormat)),ng,c(cormat))
    cordata <- subset(cordata,cordata[,1]>=cordata[,2])

    ## write to disk 
    write.table(cordata, file = paste(outheader, ".grm", sep=""), col.names=F, row.names=F, quote=F)
    system(paste("gzip ", outheader, ".grm ", sep=""))

    ## write id file to disk
    write.table(cordata.id, file = paste(outheader, ".grm.id", sep=""), col.names=F, row.names=F, quote=F)
    return(cormat)
}
