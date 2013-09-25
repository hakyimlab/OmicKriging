##
writecormat <- function(cormat,nitem,outdir,outheader){
    '%&%' <- function(a, b) paste(a, b, sep="")
    cordata <- cbind(c(col(cormat)),c(row(cormat)),nitem,c(cormat))
    cordata <- subset(cordata,cordata[,1]>=cordata[,2])

    ## write to disk 
    write.table(cordata, file = outdir %&% outheader %&% '.grm', col.names=F, row.names=F, quote=F)
    system('gzip '%&% outdir %&% outheader %&%'.grm')
}
