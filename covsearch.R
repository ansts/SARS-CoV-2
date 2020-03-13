require(Biostrings)
require(seqinr)
require(stringi)
require(stringdist)

p=224087/1.28e9

CoVprots=read.fasta("GCF_009858895.2_ASM985889v3_protein.faa")
length(CoVprots)
CovP=sapply(CoVprots,function(p){
  p=toupper(p)
  p=sapply(1:(length(p)-6), function(j){
    paste(p[j:(j+6)], collapse = "")
  })
  return(p)
})
CoVpep=unlist(CovP)
mims=read.table(file="mims.txt", header = F, quote = NULL, stringsAsFactors = F)
mimhits=mims[mims[,1] %in% CoVpep,1]
loc=sapply(CovP,function(p){
  which(p %in% mimhits)
})
loc=loc[lengths(loc)!=0]
N0=length(CoVpep)
binom.test(length(mimhits),N0,p)

spike=toupper(paste(CoVprots$YP_009724390.1, collapse=""))
x=CoVprots$YP_009724389.1[3500:3550]
x=toupper(paste(x, collapse=""))

########## SARS

SARSprots=read.fasta("SARS.fasta")
length(SARSprots)
SARSP=sapply(SARSprots,function(p){
  p=toupper(p)
  p=sapply(1:(length(p)-6), function(j){
    paste(p[j:(j+6)], collapse = "")
  })
  return(p)
})
SARSpep=unlist(SARSP)
mimSARShits=mims[mims[,1] %in% SARSpep,1]
loc1=sapply(SARSP,function(p){
  which(p %in% mimSARShits)
})
loc1=loc1[lengths(loc1)!=0]
N1=length(SARSpep)
binom.test(length(mimSARShits),N1,p)

spikepep=SARSP[[3]]

########## 229E

E229prots=read.fasta("HCoV229E.fasta")
length(E229prots)
E229P=sapply(E229prots,function(p){
  p=toupper(p)
  p=sapply(1:(length(p)-6), function(j){
    paste(p[j:(j+6)], collapse = "")
  })
  return(p)
})
E229pep=unlist(E229P)
mim229Ehits=mims[mims[,1] %in% E229pep,1]
loc2=sapply(E229P,function(p){
  which(p %in% mim229Ehits)
})
loc2=loc2[lengths(loc2)!=0]
N2=length(E229pep)
binom.test(length(mim229Ehits),N2,p)

########## OC43

OC43prots=read.fasta("OC43.fasta")
length(OC43prots)
OC43P=sapply(OC43prots,function(p){
  p=toupper(p)
  p=sapply(1:(length(p)-6), function(j){
    paste(p[j:(j+6)], collapse = "")
  })
  return(p)
})
OC43pep=unlist(OC43P)
mimOC43hits=mims[mims[,1] %in% OC43pep,1]
loc3=sapply(OC43P,function(p){
  which(p %in% mimOC43hits)
})
loc3=loc3[lengths(loc3)!=0]
N3=length(OC43pep)
binom.test(length(mimOC43hits),N3,p)

########## Poliio Mahoney

Mahprot=read.fasta("PolioMahoney.fasta")
Mahprot=toupper(Mahprot[[1]])
MahP=sapply(1:(length(Mahprot)-6), function(j){
    paste(Mahprot[j:(j+6)], collapse = "")
  })

#mimMahits=mims[mims[,1] %in% MahP,1]
mimMahits=sapply(MahP, function(p){
  p0=stringdist(mims[,1], p, method = "hamming")
  return(mims[p0<2,1])
})

mimMahits=mimMahits[lengths(mimMahits)>0]
#loc3=sapply(OC43P,function(p){
#  which(p %in% mimOC43hits)
#})
#loc3=loc3[lengths(loc3)!=0]
#N3=length(OC43pep)
#binom.test(length(mimOC43hits),N3,p)


###### Compose a synopsis table
allhits=c(mimhits,mimSARShits,mim229Ehits,mimOC43hits)
allhits=unique(allhits)
covmimpeps=unique(c(mimhits,mimSARShits,mim229Ehits,mimOC43hits))
covmpi=sapply(covmimpeps,function(p){which(allhits==p)})
locs=c(loc,loc1,loc2,loc3)
plocs=names(locs)
covproteomes=c(CoVprots,SARSprots,E229prots,OC43prots)
covpeps=c(CovP,SARSP,E229P,OC43P)

annots=sapply(plocs,function(p){attributes(covproteomes[p][[1]])$Annot})
ani=stri_extract_all(annots, regex="(?<=\\[protein=)\\w+")
anj=stri_extract_all(annots, regex="(?<=\\[protein_id=)\\w+")
anji=paste(ani,anj, sep="_")
an3=stri_extract_all(annots[1:3], regex="(?<=\\>)(.*)(?=\\[)")
anji[1:3]=an3
anji=rep(anji, times=lengths(locs))
anvir=c(rep("SARSCoV2",3),rep("SARS",4), rep('229E',4),rep("OC43",1))
anvir=rep(anvir,times=lengths(locs))
allh1=sapply(seq_along(locs), function(l){
  covpeps[names(locs)[l]][[1]][locs[[l]]]
})

CoVmimhitsdf=data.frame(Mimotopes=unlist(allh1), Protein_ID=unlist(anji), Protein=0, Strain=unlist(anvir), Starting_pos=0, stringsAsFactors = F)
CoVmimhitsdf$Starting_pos=c(unlist(loc), unlist(loc1),unlist(loc2),unlist(loc3))
CoVmimhitsdf$Protein=c("orf1ab","orf1ab","spike","orf1a","replicase","replicase","orf1ab","orf1ab","replicase","replicase","orf1a","spike","spike","hypothetical", "replicase","replicase","spike","M protein", "replicase")
write.csv(CoVmimhitsdf, file="covmimhitsdf.csv")
CoVmimhitsdf0=read.csv(file="covmimhitsdf.csv")
CoVmimhitsdf=CoVmimhitsdf0[,-1]
View(CoVmimhitsdf)

CoVmimhitsdfs=CoVmimhitsdf[-c(4,9,10,11,16),]
rownames(CoVmimhitsdfs)=1:14

IEDBSARSpike=read.csv("epitope_table_export_1583310760.csv")
IEDBepis=unique(IEDBSARSpike$Description)

IEDBlookup=lapply(CoVmimhitsdfs$Mimotopes, function(p) {
  x=stri_locate_all(IEDBepis, fixed = p)  
  x=matrix(unlist(x), byrow = T, ncol=2)
  x=x[!is.na(x[,1]),]
  return(x)
})
IEDBSARSpike[which(!is.na(stri_locate(IEDBSARSpike$Description, fixed=CoVmimhitsdfs$Mimotopes[8])[,1])),]
      

CoV2Bepipred=read.csv("SArSCoV2Bepipred2.csv")
SARSBepipred=read.csv("SArSBepipred2.csv")
CoV2Bepipred$Position=CoV2Bepipred$Position +1
SARSBepipred$Position=SARSBepipred$Position +1
SARSp_b_e=cbind(IEDBSARSpike$Starting.Position,IEDBSARSpike$Ending.Position)
SARSp_b_eo=SARSp_b_e[order(SARSp_b_e[,1]),]
plot(SARSBepipred$Score, cex=0)
lines(SARSBepipred$Score)
for (i in seq_along(SARSp_b_e[,1])){
  lines(SARSp_b_e[i,1]:SARSp_b_e[i,2],SARSBepipred$Score[SARSp_b_e[i,1]:SARSp_b_e[i,2]], col=2)
}

plot(c(0,max(SARSp_b_eo[,2])),c(0,length(SARSp_b_eo[,1])))
for (i in seq_along(SARSp_b_eo[,1])){
    lines(c(SARSp_b_eo[i,1],SARSp_b_eo[i,2]),c(i,i), col=2)
}



SARSBepipat=t(apply(SARSpatser,1,function(l){
  x=SARSBepipred$Score[SARSBepipred$Position==l[1]]
  return(c(l,x))
}))
colnames(SARSBepipat)=c("Position","SeraPercent","BepiScore")
SARSBepipat=unique(SARSBepipat)
ic=rep(1,nrow(SARSBepipat))
ic[c(53,124)]=2
labs=CoVmimhitsdfs$Mimotopes[8:9]
labs[1]=paste(389,labs[1],395, sep=" ", collapse="")
labs[2]=paste(922,labs[2],928, sep=" ", collapse="")
plot(SARSBepipat[,c(3:2)], pch=16, col=ic, xlab="Score Bepipred", ylab="% Reactive Patients' Sera ")
text(SARSBepipat[ic==2,c(3:2)], pos=2, labels=labs, col=2)
TTLDSKTrelSES=read.table("TTLDSKTrelSES.txt", header = F)
TTLDSKTrelSES=data.frame(aa=c("T","T","L","D","S","K","T"),relSES=TTLDSKTrelSES$V2)
barplot(TTLDSKTrelSES$relSES*100, names=TTLDSKTrelSES$aa, ylim=c(0,100), ylab = "% Solvent Exposed Surface")
par(new=T)
lines(c(0,10),c(5,5))

