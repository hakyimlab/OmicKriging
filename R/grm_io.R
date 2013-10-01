#' Read the GRM binary file
#'
#'
#'
#'
#'
#'
#'
#'
ReadGRMBin <- function(prefix, AllN = FALSE, size = 4, id = NULL){
  sum_i <- function(i){
    return(sum(1:i))
  }
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb");
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile <- file(NFileName, "rb");
  if(AllN==T){
    N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N <- readBin(NFile, n=1, what=numeric(0), size=size)
  i <- sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

#' Write GRM binary files
#'
#'
#'
#'
#'
#'
writeGRMBin <- function {



}

