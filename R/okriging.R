okriging <- function(idtest,idtrain=NULL,corlist,H2vec,pheno,phenoname,Xcova=NULL){
  ## corlist contains list of correlation matrices
  ## rownames and colnames of cor should be IID list and include idtest and idtrain
  ## grm is genetic relatedness matrix
  ## xrm is the expression relatedness matrix
  ## H2vec has weights for each RM relatednes matrix
  ## sig2 is total variability of phenotype
  ##
  ## universal kriging formula
  ## lambda' = ( c + X m )' iSig
  ## m' = ( x - X' iSig c )' ( X' iSig X )^-1
  ## m' = ( t(x) - c' iSig X ) ( X' iSig X )^-1
  ## lambda' = (c' + m' X) iSig
  ## x: #covariates x ntest
  ## X: ntrain x #cov
  ## c: ntrain x ntest
  ## Cressie 1993 Statistics for Spatial Data p.154

  idtest <- as.character(idtest)
  idtrain <- as.character(idtrain)
  nt <- length(idtest)
  nT <- length(idtrain)
  indall <- c(idtrain,idtest)
  if(length(unique(idtest))!=nt) warning('repeated test ids')
  if(length(unique(idtrain))!=nT) warning('repeated train ids')
  if(length(intersect(idtest,idtrain)>0)) warning('test id in training set')
  if(sum(H2vec<0) | sum(H2vec)>1) stop(' sum of weights > 1 or negative weights ')
  
  ## compute correlation matrix
  if(length(corlist)!=length(H2vec)) stop('number of correlation components (length(H2vec)) != number of corlist ')
  id <- diag(rep(1,nt+nT)) ## identity matrix
  Sigmall <- id * (1 - sum(H2vec))
  for(cc in 1:length(corlist)) Sigmall = Sigmall + H2vec[cc] * corlist[[cc]][indall,indall]

  ## row and colnames of cor should be IID
  if(sum(c(idtest,idtrain) %in% rownames(Sigmall))<(nt+nT)) stop('some correlations are missing')
  
  ## if no covariates, use intercept
  Xtest <- matrix(1,1,nt)
  Xtrain <- matrix(1,nT,1)

  if(!is.null(Xcova)) 
  {
    Xtest <- rbind(Xtest,matrix(t(Xcova[idtest,]),ncol(Xcova),nt))
    Xtrain <- cbind(Xtrain,as.matrix(Xcova[idtrain,]))
  }
  
  Ytrain <- pheno[idtrain,phenoname]
  
  ## iSig
  iSig <- solve( Sigmall[idtrain,idtrain] ) 
  
  ctvec <- matrix(Sigmall[idtest,idtrain],nt,nT )  ## correlation between new id and old id (nT x nt)
  cvec <- t(ctvec)
  mtvec <- (t(Xtest) - ctvec %*% iSig %*% Xtrain) %*% solve( t(Xtrain) %*% iSig %*% Xtrain)

  lambt <- (ctvec + mtvec %*% t(Xtrain)) %*% iSig
  Ypred <- lambt %*% Ytrain
  Ytest <- pheno[idtest,phenoname]
  res <- data.frame(IID=idtest,Ypred,Ytest)
  rownames(res) <- idtest
  return(res)
}
