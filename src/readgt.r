##
readgt <- function(genotypefile){
    ## scan a genotype matrix
    "%&%" <- function(a, b) paste(a, b, sep="")
    snpnum <- scan(pipe("awk '{print NF;exit}' "%&% genotypefile))
    idnum <- scan(pipe("grep -c ^ "%&% genotypefile)) - 1
    gtmat <- matrix(scan(genotypefile, n = snpnum * idnum,
    skip = 1, na.strings = "NA", what = "character"), idnum, snpnum,
    byrow=TRUE)
    print('Read genotype matrix')
    dims <-t dim(gtmat)
    gtmat <- gtmat[,7:dims[2]] 
    gtmat <- as.numeric(gtmat)
    nas <- which(is.na(gtmat), arr.ind=TRUE)
    gtmat[nas] <- colMeans(gtmat, na.rm=TRUE)[nas[,2]]
    return(gtmat)
}
