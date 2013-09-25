readcorlist <-
function(corfilelist)
{
  nc = length(corfilelist)
  corlist = list()
  for(cc in 1:nc)
  {
    corlist[[cc]] = grm2mat(corfilelist[cc])
  }
  return(corlist)
}
