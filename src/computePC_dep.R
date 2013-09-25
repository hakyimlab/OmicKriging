computePC <-
function(fullheader,idfile=idfile,gctaname="gcta")
{
  runPC = paste(gctaname , ' --grm ' , fullheader , ' --keep ' , idfile , ' --pca 10 --out ' , fullheader,sep="")
  print(runPC)
  system(runPC)
}
