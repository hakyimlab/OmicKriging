##
readpheno <- function(phenofile){
    pheno <- read.table(phenofile,as.is=T,header=F)
    names(pheno) <- c("FID","IID","Status")
    rownames(pheno) <- pheno$IID
    return(pheno)
}
