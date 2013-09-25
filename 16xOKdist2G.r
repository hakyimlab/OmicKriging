## ============================
## OK
## ============================

args <- commandArgs(trailingOnly=TRUE)

library(ROCR)
library(parallel)
library(doMC)
library(OmicKriging)
library(MASS)

okrigingfast <- function (idtest, idtrain = NULL, corlist, H2vec, pheno, phenoname, 
    Xcova = NULL) 
{
    idtest = as.character(idtest)
    idtrain = as.character(idtrain)
    nt = length(idtest)
    nT = length(idtrain)
    indall = c(idtrain, idtest)
    if (length(unique(idtest)) != nt) 
        warning("repeated test ids")
    if (length(unique(idtrain)) != nT) 
        warning("repeated train ids")
    if (length(intersect(idtest, idtrain) > 0)) 
        warning("test id in training set")
    if (sum(H2vec < 0) | sum(H2vec) > 1) 
        stop(" sum of weights > 1 or negative weights ")
    if (length(corlist) != length(H2vec)) 
        stop("number of correlation components (length(H2vec)) != number of corlist ")
    id = diag(rep(1, nt + nT))
    Sigmall = id * (1 - sum(H2vec))
    for (cc in 1:length(corlist)) Sigmall = Sigmall + H2vec[cc] * 
        corlist[[cc]][indall, indall]
    if (sum(c(idtest, idtrain) %in% rownames(Sigmall)) < (nt + 
        nT)) 
        stop("some correlations are missing")
    Xtest = matrix(1, 1, nt)
    Xtrain = matrix(1, nT, 1)
    if (!is.null(Xcova)) {
        Xtest = rbind(Xtest, matrix(t(Xcova[idtest, ]), ncol(Xcova), 
            nt))
        Xtrain = cbind(Xtrain, as.matrix(Xcova[idtrain, ]))
    }
    Ytrain = pheno[idtrain, phenoname]
    iSig = chol2inv(chol(Sigmall[idtrain, idtrain]))
    ctvec = matrix(Sigmall[idtest, idtrain], nt, nT)
    cvec = t(ctvec)
    newvec = t(Xtrain) %*% iSig %*% Xtrain
    mtvec = (t(Xtest) - ctvec %*% iSig %*% Xtrain) %*% chol2inv(chol(newvec))
    lambt = (ctvec + mtvec %*% t(Xtrain)) %*% iSig
    Ypred = lambt %*% Ytrain
    Ytest = pheno[idtest, phenoname]
    res = data.frame(IID = idtest, Ypred, Ytest)
    rownames(res) = idtest
    return(res)
}

foldername <- args[1]
ndist <- 1:args[2]
H2vec <- c(as.numeric(args[3]),as.numeric(args[4]))

"%&%" <- function(a, b) paste(a, b, sep="")

pre <- '/nas40t0/keston/OKt2/'
tempo.dir <- pre %&% 'tempo/' 

work.dir = "."
setwd(work.dir)

## files
phenofile <- tempo.dir %&% foldername %&% "Pheno"
grmfullheader <- tempo.dir %&% foldername
grmpartheader <- tempo.dir %&% foldername %&% "_S"


## read in grm
cortempo <- data.frame(headers=c(grmfullheader,grmpartheader),stringsAsFactors=F)
corfilelist <- c(cortempo[,1])
corlist <- readcorlist(corfilelist) 

## read pheno
pheno <- read.table(phenofile,as.is=T,header=F)
names(pheno) <- c("FID","IID","Status")
rownames(pheno) <- pheno$IID
phenoname <- "Status"

phenovec <- as.vector(as.numeric((pheno[,3])))
idvec <- 1:length(phenovec)
names(phenovec) <- idvec

iddata <- pheno[,1:3]
names(iddata) <- c("FID","IID","Status")

# setup the cores
ncore <- detectCores()
registerDoMC(cores = ncore)

#output <- data.frame(rep(NA, length(ndist)),ndist,foldername)
#names(output) <- c("AUC","Dist","WT")

for(x in ndist){
    print('Running -> '%&% x)
    ##
    idlength <- length(iddata$IID)
    groupid <- sample(1:ncore, idlength, replace=T)

    newiddata <- data.frame(groupid,iddata$IID)
    colnames(newiddata) <- c("GID","IID")
    idsubsetlist <- 1:ncore

    # run kriging in parallel
    p <- foreach(i=1:ncore,.combine=rbind) %dopar% {
    idtest <- newiddata$IID[newiddata$GID == idsubsetlist[i]]
    idtrain <- newiddata$IID[!(newiddata$IID %in% idtest)]
    okriging(idtest,idtrain,corlist,H2vec,pheno,phenoname=phenoname)
    }
    
    gc()

    # summary
    pred <- prediction(p$Ypred, p$Ytest)
    perf <- performance(pred,"auc")
    print(perf@y.values)
    #output[x,1] <- perf@y.values
    write.table(perf@y.values,file = pre %&% foldername %&% '_finalv2_appendout',quote=F,col.names=F,row.names=F,append=T)
}
#write.table(output,file = pre %&% foldername %&% '_2gx_output',quote=F,col.names=T,row.names=F)
## END
