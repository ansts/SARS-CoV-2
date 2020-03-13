# Mimotopes

hALmm=pepD5(mims[,1])
Ghamm=adjL2G(hALmm)
vgc=components(Ghamm)$membership==1
save(Ghamm, file="Ghamm_graph")
Ghmm=induced.subgraph(Ghamm, vgc)
rm(Ghamm)

# adding weights
vGhmm=attributes(V(Ghmm))$names
eGhmm=attributes(E(Ghmm))$vnames
whALmm=sapply(eGhmm,function(n){
  l=unlist(strsplit(n, split="\\|"))
  x=stringdist(l[1],l[2], method="hamming")
  return(x)
})

edge.attributes(Ghmm)$weight=1/whALmm

write.graph(Ghmm, file = "Ghmmw.graphml", format = "graphml")

covmimtable=read.csv("covmimtable.csv")
mimhitsall=unique(covmimtable$Mimotopes)
miia=vGhmm[vGhmm %in% mimhitsall]
mimhitsallego=make_ego_graph(Ghmm, nodes=miia, mode="all")
mimhitsegos=lapply(mimhitsallego, function(g){
  attributes(V(g))$names
})

BinT=sapply(lengths(mimhitsegos), function(n){
  binom.test(n, 19*7+19^2*choose(7,2), p, alternative = "g")$p.value
})
BinTa=cbind(lengths(mimhitsegos), p.adjust(BinT, method="BH"))
