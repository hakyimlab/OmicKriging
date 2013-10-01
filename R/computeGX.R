#' Compute gene expression correlation matrix.
#'
#'
#'
#'
#'
#'
#'
#'
computeGX <- function(genefile,outheader,idfile=NULL) {
  genedata = read.delim(genefile,sep="",as.is=T,header=T)
  if(!is.null(idfile)) {
    iddata = read.table(idfile,header=T,as.is=T)
    genedata = merge(iddata,genedata,by.x=c("FID","IID"),by.y=c("FID","IID"))
  }
  cordata.id = genedata[c("FID","IID")]
  genemat = as.matrix(genedata[,!(names(genedata) %in% c("FID","IID"))])
  ## impute missing data
  ## Standardize this matrix?
  if(sum(is.na(genemat))>0) print("Impute first") ## TODO: impute with the script
  cormat = cor(t(genemat))
  ## format matrix into linear dataframe (GCTA format)
  ng = ncol(genemat)
  cordata = cbind(c(col(cormat)),c(row(cormat)),ng,c(cormat))
  cordata = subset(cordata,cordata[,1]>=cordata[,2])
  ## write to disk
  ## replace this with write to GCTA style correlation matrix file
  write.table(cordata,file=paste(outheader , ".grm",sep=""),col.names=F,row.names=F,quote=F)
  system(paste("gzip " , outheader , ".grm ",sep=""))
  ## write id file to disk
  write.table(cordata.id,file=paste(outheader , ".grm.id",sep=""),col.names=F,row.names=F,quote=F)
  return(cormat)
}
