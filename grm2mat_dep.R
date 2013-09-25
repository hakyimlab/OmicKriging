grm2mat <-
function(grmfullheader)
{
  grm = read.table(paste(grmfullheader , ".grm.gz",sep=""),header=F,as.is=T)
  grmid = read.table(paste(grmfullheader , ".grm.id",sep=""),header=F,as.is=T)
  names(grmid) = c("FID","IID")
  grmid$IID = as.character(grmid$IID)
  nn = nrow(grmid)
  grmat = matrix(0,nn,nn)
  ## computes linear index of matrix corresponding to ii,jj
  indexfun = function(i,j,nn) return(i+(j-1)*nn)
  grmat[indexfun(grm[,1],grm[,2],nn) ] = grm[,4]  
  grmat= t( lower.tri(grmat) * grmat ) + grmat  ## 
  rownames(grmat) = grmid$IID
  colnames(grmat) = grmid$IID
  return(grmat)
}
