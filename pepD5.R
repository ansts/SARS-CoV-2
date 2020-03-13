
pepD5=function(p){
  require(Biostrings)
  require(parallel)

  cl = makeCluster(4)
  clusterExport(cl,varlist=c("p"), envir = environment())
  clusterEvalQ(cl, require(stringdist))
  
  d=parSapply(cl,seq_along(p),function(i){
       pi=p[i]
       x=stringdist(p[-i],pi,method = "hamming")
       m=p[-i][x<3]
    return(m)
  })
  stopCluster(cl)
  names(d)=p
  return(d)
}