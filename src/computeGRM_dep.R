computeGRM <-
function(bedfullheader,grmfullheader="genomic",gctaname="gcta",idfile=NULL,snpfile=NULL)
{
  if(is.null(idfile)) idfile = paste(bedfullheader , ".fam",sep="")
  runGRM = paste(gctaname , ' --bfile ' , bedfullheader , ' --autosome --make-grm --keep ' , idfile , sep="") 
  if(!is.null(snpfile)) runGRM = paste(runGRM , '--extract ' , snpfile , sep="")
  runGRM = paste(runGRM , ' --out ' , grmfullheader, sep="")
  system(runGRM)
}
