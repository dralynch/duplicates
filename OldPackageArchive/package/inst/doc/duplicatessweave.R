### R code from vignette source 'duplicatessweave.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()
library("LynchSmithEldridgeTavareFragDup")
endloop<-22


###################################################
### code chunk number 2: Overview
###################################################
list.files(system.file("extdata", package="LynchSmithEldridgeTavareFragDup"))


###################################################
### code chunk number 3: Preparation
###################################################
library("Rsamtools")
library("BSgenome.Hsapiens.UCSC.hg19")
EDpath <- system.file("extdata", package="LynchSmithEldridgeTavareFragDup")


###################################################
### code chunk number 4: PicardUnmasked
###################################################
unmaskedP<-read.delim(file.path(EDpath,"Picard", "samples.metrics.txt"),as.is=T)

readPairsExamined <- 
  readPairDuplicates <- 
  OpticalPairDuplicates <- matrix(ncol=9, nrow=22, 
                            dimnames = list(unmaskedP[,1], c("Total", "X", "Y", "M",
                                                             "Centromeres","Telomeres",
                                                             "LowCov","HighCov",
                                                             "Residual")))

readPairsExamined[,1]<-unmaskedP[,3]
readPairDuplicates[,1]<-unmaskedP[,6]
OpticalPairDuplicates[,1]<-unmaskedP[,7]


###################################################
### code chunk number 5: LatexTableValues
###################################################
latexdiv22<-rep("&",22)
latexend22<-rep("\\\\",22)
latexdiv8<-rep("&",8)
latexend8<-rep("\\\\",8)


###################################################
### code chunk number 6: ReadWeaver
###################################################
WeaverSuppTable1<-read.delim(file.path(EDpath,"WeaverSuppTable1.txt"),sep="",header=T,as.is=T)
WeaverSuppTable1$CaseID<-WeaverSuppTable1$CaseID-1


###################################################
### code chunk number 7: GenTable1
###################################################
Table1<-cbind(
WeaverSuppTable1[,1],latexdiv22,
WeaverSuppTable1[,2],latexdiv22,
WeaverSuppTable1[,12],latexdiv22,
round(unmaskedP[,3]/1000000000,2),latexdiv22,
round(100*unmaskedP[,6]/unmaskedP[,3],2),latexdiv22,
round(unmaskedP[,9]/1000000000,2),latexend22
)

Table1[Table1=="NormalOesophagus"]<-"NormalOes"

Table1<-rbind(c("ID", "&","Sex","&","Tissue","&","Reads","&","Dup.","&","Library", "\\\\"),
              c("","&","","&","Tissue","&","($10^9$)","&","Rate","&","Size","\\\\" ), 
              c("","&","","&","","&","","&","","&","($10^9$)","\\\\"), Table1)

write.table(Table1,sep=" ",file="Table1.tsv",row.names=F,col.names=F,quote=F)
Table1


###################################################
### code chunk number 8: ExploreSI2
###################################################
table(WeaverSuppTable1[,2],WeaverSuppTable1[,12])
Group<-1+2*(WeaverSuppTable1[,2]=="Male")+(WeaverSuppTable1[,12]=="Blood")
GroupN<-c("Female_Tissue","Female_Blood","Male_Tissue","Male_Blood")[Group]


###################################################
### code chunk number 9: exploreMasks
###################################################
list.files(file.path(EDpath,"masks"))
head(read.delim(file.path(EDpath, "masks", "Mask4-Centromeres.bed"),header=F,as.is=T))

maskFiles <- list.files(file.path(EDpath,"masks"), full.names = TRUE)
MASKS <- lapply(maskFiles, read.delim, header=FALSE, as.is=TRUE, skip=1)
effectiveGenomeSize<-rep(0,9)
effectiveGenomeSize[1]<-sum(as.numeric(MASKS[[5]][17:38,3]))
+MASKS[[1]][1,3]+MASKS[[2]][1,3]+MASKS[[3]][1,3]
for(i in 1:7){ 
    effectiveGenomeSize[i+1] <- sum(MASKS[[i]][,3] - MASKS[[i]][,2])}
effectiveGenomeSize[9]<-effectiveGenomeSize[1]-sum(effectiveGenomeSize[2:8])


###################################################
### code chunk number 10: PicardMasked
###################################################
metfiles<-list.files(file.path(EDpath,"Picard"))
metfiles<-grep("mask",metfiles,value=T)

for(i in 1:7){
  temp<-read.delim(file.path(EDpath,"Picard",metfiles[i]),as.is=T)
  temp<-temp[-(15:19),]
  readPairsExamined[,i+1]<-temp[,3]
  readPairDuplicates[,i+1]<-temp[,6]
  OpticalPairDuplicates[,i+1]<-temp[,7]
}


###################################################
### code chunk number 11: PicardResid
###################################################
readPairsExamined[,9]<-readPairsExamined[,1]-apply(readPairsExamined[,2:8],1,sum)
readPairDuplicates[,9]<-readPairDuplicates[,1]-apply(readPairDuplicates[,2:8],1,sum)
OpticalPairDuplicates[,9]<-OpticalPairDuplicates[,1]-apply(OpticalPairDuplicates[,2:8],1,sum)


###################################################
### code chunk number 12: PicardTable2
###################################################

reptable<-matrix(0,ncol=4,nrow=8)
reptable[1,]<-round(100*sapply(split(((readPairDuplicates[,9]/readPairsExamined[,9])),Group),mean),2)
for(i in 2:8){
  reptable[i,]<-round(sapply(split((
    (readPairDuplicates[,i]/readPairsExamined[,i])
    /(readPairDuplicates[,9]/readPairsExamined[,9])),Group),mean),2)
}

Table2<-cbind(reptable[,1],latexdiv8,reptable[,2],latexdiv8,reptable[,3],latexdiv8,
              reptable[,4],latexend8)
write.table(Table2,sep="\t",file="Table2.tsv",row.names=F,col.names=F,quote=F)

reptable[3,1:2]<-"-"
rownames(reptable)<-c("Residual","X","Y","M","Centromeres","Telomeres","Low Cov","High Cov")
colnames(reptable)<-c("Female Tissue","Female Blood","Male Tissue","Male Blood")
reptable


###################################################
### code chunk number 13: LoadSNPlist
###################################################
snplist<-read.delim(file.path(EDpath,"SNPstoextract.txt"),as.is=T,header=T)


###################################################
### code chunk number 14: AddFig1
###################################################
par(mar=c(4.6,4.1,1.6,1.6))
plot(table(snplist$chrom)[c(1,12,16:22,2:11,13:15,23:24)], 
     ylab="number of SNPs",xlab="chromosome",axes=F)
axis(2)
axis(1,labels=c(1:22,"X","Y"),at=1:24,las=2)
box()


###################################################
### code chunk number 15: SNPGCbias
###################################################
genome <- Hsapiens
GCNo<-rep(NA,2500)
for(i in 1:2500){
  GCNo[i]<-sum(unlist(strsplit(as.character(substr(genome[[snplist$chrom[i]]],snplist$chromEnd[i]-250,snplist$chromEnd[i]+249)),"")) %in% c("G","C"))
}


###################################################
### code chunk number 16: GCrefgenome
###################################################
GCbinned<-NULL
for(k in 1:23){
    windowViews <- trim(Views(genome[[k]], start = seq(1, length(genome[[k]]), 500), width = 500))
    letterFreq <- letterFrequency(windowViews, letters = c("A","C","G","T"))
    seqGC <- rowSums(letterFreq[,2:3]) / rowSums(letterFreq)
    GCbinned<-c(GCbinned,seqGC)
}


###################################################
### code chunk number 17: AddFigSNPGC
###################################################
par(mfrow=c(1,2))
hist(GCbinned,freq=F,col=rgb(1,0.5,0.5,0.5),breaks=seq(0,1,0.02),main="Distribution of GC content of genome",xlab="GC proportion")
hist(GCNo/500,freq=F,col=rgb(0.5,0.5,1,0.5),breaks=seq(0,1,0.02),add=T)
legend("topright",fill=c(rgb(1,0.5,0.5,0.5),rgb(0.5,0.5,1,0.5)),legend=c("Genome-wide","Our SNPs"))

plot(quantile(GCbinned,probs=seq(0.01,.99,0.01),na.rm=T),quantile(GCNo/500,probs=seq(.01,.99,0.01),na.rm=T),
     pch=20,col="red",xlab="Genome wide 1st - 99th percentiles",
     ylab="Our SNPs 1st - 99th percentiles",
     main="Proportion GC Content\nKolmogorov-Smirnov p=0.08")
for(i in 1:99){
  GW<-quantile(GCbinned,probs=i/100,na.rm=T)
  OS<-quantile(GCNo/500,probs=i/100,na.rm=T)
lines(c(GW,GW),c(0,OS),lwd=0.5)  
lines(c(0,GW),c(OS,OS),lwd=0.5)  
}
points(quantile(GCbinned,probs=seq(0.01,.99,0.01),na.rm=T),quantile(GCNo/500,probs=seq(.01,.99,0.01),na.rm=T),pch=20,col="red")
abline(0,1,lwd=2)


###################################################
### code chunk number 18: KStest
###################################################
ks.test(GCbinned,GCNo/500)


###################################################
### code chunk number 19: Het SNP processing1 (eval = FALSE)
###################################################
## bamfilelist <- list.files("data/snpbams", pattern = ".bam$")


###################################################
### code chunk number 20: Het SNP processing2 (eval = FALSE)
###################################################
## ACvec<-c("2:0", "1:1" ,
##          "3:0" ,"2:1" ,
##          "4:0" ,"3:1" ,"2:2" ,
##          "5:0" ,"4:1" ,"3:2" ,
##          "6:0" ,"5:1" ,"4:2" ,"3:3" ,
##          "7:0" ,"6:1" ,"5:2" ,"4:3")


###################################################
### code chunk number 21: Het SNP processing3 (eval = FALSE)
###################################################
## HetSNPTable<-matrix(0,nrow=length(bamfilelist),ncol=18)
## SNPnumbers<-rep(0,22)
## HSTrow<-0


###################################################
### code chunk number 22: Het SNP processing3b (eval = FALSE)
###################################################
## HetSNPTableHighGC<-matrix(0,nrow=length(bamfilelist),ncol=18)
## HetSNPTableLowGC<-matrix(0,nrow=length(bamfilelist),ncol=18)
## SNPnumbersHigh<-rep(0,22)
## SNPnumbersLow<-rep(0,22)


###################################################
### code chunk number 23: Het SNP processing4 (eval = FALSE)
###################################################
## for(sample in bamfilelist){
##   HSTrow<-HSTrow+1
##   cat(sample,"\n")
##   bamfile=paste("../finaldupsanalysis/data/snpbams/",sample,sep="")
##   
##   # look to see which of the candidate SNPs is heterozygous in this sample. 
##   # Note that this is based on all the reads not just the duplicate fragments, 
##   # so should not greatly bias matters.
##   
##   fls <- PileupFiles(bamfile)  
##   which<-GRanges(snplist[1:2500,1],IRanges(snplist[1:2500,3],snplist[1:2500,3]))
##   PUP <- ApplyPileupsParam(which=which, yieldSize=1000000L, yieldBy="position", what="seq",maxDepth=200,minDepth=0)
##   outres <- applyPileups(fls,(function(x){x[["seq"]]}),param=PUP) 
##   
##   # Our requirement is that the minor allele frequency is greater than 0.4 - 
##   # quite a stringent criterion  
##   
##   usesnplist<-snplist[
##     which(apply((outres)[[1]],3,secondbiggest)/apply((outres)[[1]],3,sum)>0.4),]
##   
##     useGCNo<-GCNo[which(apply((outres)[[1]],3,secondbiggest)/apply((outres)[[1]],3,sum)>0.4)]
##     SNPnumbersHigh[HSTrow]<-sum(useGCNo>=199)
##     SNPnumbersLow[HSTrow]<-sum(useGCNo<199)
## 
##   
##   SNPnumbers[HSTrow]<-dim(usesnplist)[1]
##   # now we are going to count up the numbers of duplicates that share the 
##   # same allele and the numbers that do not.
##   
##   truetally<-0
##   falsetally<-0
##   
##   store<-rep(0,dim(usesnplist)[1])
##   # for each allele on the list
##   for(snp in 1:dim(usesnplist)[1]){
##     newwhich<-GRanges(usesnplist[snp,1], IRanges(usesnplist[snp,3], usesnplist[snp,3]))
##     
##     #these are the reads that were marked as duplicates 
##     ydfile<-scanBam(bamfile, 
##                     param=ScanBamParam(flag=scanBamFlag(isDuplicate=T),
##                                        simpleCigar=T,what=c("pos","mpos","seq"),
##                                        which=newwhich))[[1]]
##     
##     #these are the reads that were not
##     ndfile<-scanBam(bamfile,
##                     param=ScanBamParam(flag=scanBamFlag(isDuplicate=F),
##                                        simpleCigar=T,what=c("pos","mpos","seq"),
##                                        which=newwhich))[[1]]
##     
##     # Any fragement marked as a duplicate must be a duplicate of a fragment 
##     # that is not marked as a duplicate. We just want to keep a matched set 
##     # of duplicates and the fragments of which they are duplicates
##     
##     ykey<-paste(ydfile$pos,ydfile$mpos)
##     nkey<-paste(ndfile$pos,ndfile$mpos)
##     
##     usekey<-unique(ykey)
##     
##     # Now, assuming that we see some duplicates, we are going to go through them 
##     # and compare the alleles
##     if(length(usekey)>0){
##       for(key in usekey){
##         
##         storealleles<-NULL
##         for(p in which(ykey==key)){
##           storealleles<-c(storealleles, substr(ydfile$seq[p], 
##                                                usesnplist[snp,3]-ydfile$pos[p]+1,
##                                                usesnplist[snp,3]-ydfile$pos[p]+1))
##         }
##         for(p in which(nkey==key)){
##           storealleles<-c(storealleles,substr(ndfile$seq[p],
##                                               usesnplist[snp,3]-ndfile$pos[p]+1,
##                                               usesnplist[snp,3]-ndfile$pos[p]+1))
##         }
##         
##         ortab<-outres[[1]][,,as.numeric(rownames(usesnplist)[snp])]
##         usebase<-names(sort(ortab,decreasing=T))[1:2]
##         
##         
##         allelecounts<-c(sum(storealleles==usebase[1]),sum(storealleles==usebase[2]))
##         if(allelecounts[2]>allelecounts[1]){allelecounts<-allelecounts[2:1]}
##         
##         HetSNPTable[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]<-
##           HetSNPTable[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]+1
##         if(useGCNo[snp]>=199){
##           HetSNPTableHighGC[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]<-
##             HetSNPTableHighGC[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]+1}
##           
##           if(useGCNo[snp]<199){
##             HetSNPTableLowGC[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]<-
##               HetSNPTableLowGC[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]+1}
##       }
##     }
##   }
## }


###################################################
### code chunk number 24: Het SNP processing5 (eval = FALSE)
###################################################
## HetSNPTable<-cbind(SNPnumbers,HetSNPTable)
## HetSNPTableHighGC<-cbind(SNPnumbersHigh,HetSNPTableHighGC)
## HetSNPTableLowGC<-cbind(SNPnumbersLow,HetSNPTableLowGC)
## 
## colnames(HetSNPTable)<-c("NoSNPs", "AA", "AB", 
##                            "AAA", "AAB", 
##                            "AAAA", "AAAB", "AABB", 
##                            "AAAAA", "AAAAB", "AAABB", 
##                            "AAAAAA", "AAAAAB", "AAAABB", "AAABBB", 
##                            "AAAAAAA", "AAAAAAB", "AAAAABB", "AAAABBB")
## colnames(HetSNPTableLowGC)<-c("NoSNPs", "AA", "AB", 
##                          "AAA", "AAB", 
##                          "AAAA", "AAAB", "AABB", 
##                          "AAAAA", "AAAAB", "AAABB", 
##                          "AAAAAA", "AAAAAB", "AAAABB", "AAABBB", 
##                          "AAAAAAA", "AAAAAAB", "AAAAABB", "AAAABBB")
## colnames(HetSNPTableHighGC)<-c("NoSNPs", "AA", "AB", 
##                          "AAA", "AAB", 
##                          "AAAA", "AAAB", "AABB", 
##                          "AAAAA", "AAAAB", "AAABB", 
##                          "AAAAAA", "AAAAAB", "AAAABB", "AAABBB", 
##                          "AAAAAAA", "AAAAAAB", "AAAAABB", "AAAABBB")
## write.table(HetSNPTable,file="HetSNPDups.txt",sep="\t")
## write.table(HetSNPTableHighGC,file="HetSNPDupsHighGC.txt",sep="\t")
## write.table(HetSNPTableLowGC,file="HetSNPDupsLowGC.txt",sep="\t")


###################################################
### code chunk number 25: LoadSNPsummary1
###################################################
HetSNPTable<-read.delim(file.path(EDpath,"HetSNPDups.txt"),as.is=T)


###################################################
### code chunk number 26: LoadSNPsummary2
###################################################
ResDupR<-readPairDuplicates[,9]/readPairsExamined[,9]


###################################################
### code chunk number 27: LoadSNPsummary3
###################################################
ResDepth<-round(200*readPairsExamined[,9]/effectiveGenomeSize[9],2)


###################################################
### code chunk number 28: LoadSNPsummary3
###################################################
summary(HetSNPTable[,1])


###################################################
### code chunk number 29: LoadSNPsummary4
###################################################
summary(lm(HetSNPTable[,1]~ResDepth+ResDupR))
cor.test(ResDepth,HetSNPTable[,1])


###################################################
### code chunk number 30: LoadSNPsummary5
###################################################
summary(apply(HetSNPTable[,-1],1,sum))


###################################################
### code chunk number 31: LoadSNPsummary6
###################################################
summary(lm(apply(HetSNPTable[,-1],1,sum)~ResDupR+HetSNPTable[,1]+ResDepth))


###################################################
### code chunk number 32: LoadSNPsummary7
###################################################
sum(HetSNPTable[,-1])


###################################################
### code chunk number 33: LoadSNPsummary8
###################################################
sum(HetSNPTable[,2:3])
round(100*sum(HetSNPTable[,2:3])/sum(HetSNPTable[,-1]))


###################################################
### code chunk number 34: LoadSNPsummary9
###################################################
sum(HetSNPTable[,4:5])
round(100*sum(HetSNPTable[,4:5])/sum(HetSNPTable[,-1]))


###################################################
### code chunk number 35: LoadSNPsummary10
###################################################
sum(HetSNPTable[,6:8])
round(100*sum(HetSNPTable[,6:8])/sum(HetSNPTable[,-1]),1)


###################################################
### code chunk number 36: LoadSNPsummary11
###################################################
sum(HetSNPTable[,16:19])


###################################################
### code chunk number 37: GenPL
###################################################
partitionlist<-genPL(7)
CoefList<-genCoefs(partitionlist)


###################################################
### code chunk number 38: TimeFig
###################################################
PLtimes<-c(.001,.003,.007,.03,.119,.584,1.97,6.98,18.9,67.529,152.94,541.793)
CLtimes<-c(.003,.052,1.088,16.711,212.698,NA,NA,NA,NA,NA,NA,NA)
Mvals<-c(5,10,15,20,25,30,35,40,45,50,55,60)
plot(Mvals,log10(PLtimes),ylab="time (seconds)",xlab="M",axes=F)
axis(1)
axis(2,10^seq(-3,2,1),at=seq(-3,2,1),las=2)
box()
points(Mvals,log10(CLtimes),pch=17)
lm(log10(PLtimes)~Mvals)
lm(log10(CLtimes)~Mvals)

abline(-3.58,0.1071)
abline(-3.7064,0.2442,lwd=2)

legend("bottomright",pch=c(1,17),legend=c("Partition List","Coefficients"))


###################################################
### code chunk number 39: EstimateFDvec
###################################################
FDvec<-rep(0,22)
for(j in 1:endloop){
    myseq<-seq(0,0.2,length.out=1000)
    vals <- sapply(myseq, propllik, x = HetSNPTable[j,-1])
    FDvec[j]<-myseq[which.max(vals)]
}


###################################################
### code chunk number 40: EstimateDupRates
###################################################
  BasicEst<-(readPairDuplicates[,1]-OpticalPairDuplicates[,1])/
  (readPairsExamined[,1]-OpticalPairDuplicates[,1])
  ResidEst<-(readPairDuplicates[,9]-OpticalPairDuplicates[,9])/
  (readPairsExamined[,9]-OpticalPairDuplicates[,9])
  FDCorEst<-(1-FDvec)*ResidEst


###################################################
### code chunk number 41: AddFig2
###################################################
  par(mfrow=c(1,2))
  plot(c(1,3),c(0,0.15),type="n",ylab="duplicate rate",xlab="",axes=F,yaxs="i")
  axis(2)
  box()
  for(i in seq(0,.15,.03)){
    rect(-1,i-.015,5,i,border="grey95",col="grey95")
  }
  for(i in 1:22){
    tcol<-"black"
    if(WeaverSuppTable1[i,12]=="Blood"){tcol<-"red"}
    lines(1:3,c(BasicEst[i],ResidEst[i],FDCorEst[i]),type="b",pch=16,col=tcol)
  }
  axis(1,at=1:3,labels=c("Standard\n estimate","Residual\n estimate","Corrected\n estimate"),las=3)
  
  plot(ResidEst/BasicEst,FDCorEst/ResidEst, xlab="Residual/Basic estimates",
       ylab="Corrected/Residual estimates",pch=20)


###################################################
### code chunk number 42: WhenMis2
###################################################
FDvec2<-rep(0,22)

for(j in 1:endloop){
  vals<-rep(0,1000)
  myseq<-seq(0,0.2,length.out=1000)
  for(i in 1:1000){
    tmpx<-HetSNPTable[j,-1]
    tmpx[3:18]<-0
    vals[i]<-propllik(tmpx,myseq[i])
  }
  
  FDvec2[j]<-myseq[which.max(vals)]
}


###################################################
### code chunk number 43: WhenMisnot2
###################################################
FDvecN2<-rep(0,22)

for(j in 1:endloop){
  vals<-rep(0,1000)
  myseq<-seq(0,0.2,length.out=1000)
  for(i in 1:1000){
    tmpx<-HetSNPTable[j,-1]
    tmpx[1:2]<-0
    vals[i]<-propllik(tmpx,myseq[i])
  }
  
  FDvecN2[j]<-myseq[which.max(vals)]
}


###################################################
### code chunk number 44: ConcFig
###################################################
plot(FDvec2,FDvecN2,xlab="estimates from sets of 2 duplicate fragments",
     ylab="estimates from sets of >2 duplicate fragments")
abline(0,1)


###################################################
### code chunk number 45: LoadSNPsummary1
###################################################
HetSNPTableHighGC<-read.delim(file.path(EDpath,"HetSNPDupsHighGC.txt"),as.is=T)
HetSNPTableLowGC<-read.delim(file.path(EDpath,"HetSNPDupsLowGC.txt"),as.is=T)


###################################################
### code chunk number 46: GCComparison
###################################################
FDvecHigh<-rep(0,22)
FDvecLow<-rep(0,22)
for(j in 1:22){
  #cat(j,"\t")
  myseq<-seq(0,0.2,length.out=1000)
  valsHigh <- sapply(myseq, propllik, x = HetSNPTableHighGC[j,-1])
  valsLow <- sapply(myseq, propllik, x = HetSNPTableLowGC[j,-1])
  FDvecHigh[j]<-myseq[which.max(valsHigh)]
  FDvecLow[j]<-myseq[which.max(valsLow)]
}


###################################################
### code chunk number 47: GCCompFig
###################################################
  plot(FDvecHigh,FDvecLow,xlab="estimates from duplicates in 'high GC' regions",
       ylab="estimates from duplicates in 'low GC' regions")
abline(0,1)


###################################################
### code chunk number 48: LoadTumourPicardData
###################################################
TumourPicard<-read.csv(file.path(EDpath,"Tumour","Picard.csv"),as.is=T)
TumourPicard
# and the duplicate rates
round((TumourPicard[,6]-TumourPicard[,7])/(TumourPicard[,3]-TumourPicard[,7]),3)


###################################################
### code chunk number 49: TumourSNPlists
###################################################
snplistAABB<-read.delim(file.path(EDpath,"Tumour","usedupAABB.txt"),as.is=T,header=F)
snplistAABB<-snplistAABB[,c(1,2,2)]
snplistAB<-read.delim(file.path(EDpath,"Tumour","usedupAB.txt"),as.is=T,header=F)
snplistAB<-snplistAB[,c(1,2,2)]


###################################################
### code chunk number 50: ProcessTumourBam (eval = FALSE)
###################################################
## 
## HetSNPTableAABB<-matrix(0,nrow=2,ncol=18)
## bamfile="SS6003314.bam"
## usesnplist<-snplistAABB
## 
## #SNPnumbers[HSTrow]<-dim(usesnplist)[1]
## # now we are going to count up the numbers of duplicates that share the 
## # same allele and the numbers that do not.
## 
## # for each allele on the list
## for(snp in 1:dim(usesnplist)[1]){
##   if(100*trunc(snp/100)==snp){cat(snp,"\t")}
##   newwhich<-GRanges(usesnplist[snp,1], IRanges(usesnplist[snp,3], usesnplist[snp,3]))
##   
##   #these are the reads that were marked as duplicates 
##   ydfile<-scanBam(bamfile, 
##                   param=ScanBamParam(flag=scanBamFlag(isDuplicate=T),
##                                      simpleCigar=T,what=c("pos","mpos","seq"),
##                                      which=newwhich))[[1]]
##   
##   #these are the reads that were not
##   ndfile<-scanBam(bamfile,
##                   param=ScanBamParam(flag=scanBamFlag(isDuplicate=F),
##                                      simpleCigar=T,what=c("pos","mpos","seq"),
##                                      which=newwhich))[[1]]
##   
##   # Any fragement marked as a duplicate must be a duplicate of a fragment 
##   # that is not marked as a duplicate. We just want to keep a matched set 
##   # of duplicates and the fragments of which they are duplicates
##   
##   ykey<-paste(ydfile$pos,ydfile$mpos)
##   nkey<-paste(ndfile$pos,ndfile$mpos)
##   
##   usekey<-unique(ykey)
##   
##   # Now, assuming that we see some duplicates, we are going to go through them 
##   # and compare the alleles
##   if(length(usekey)>0){
##     for(key in usekey){
##       
##       storealleles<-NULL
##       for(p in which(ykey==key)){
##         storealleles<-c(storealleles, substr(ydfile$seq[p], 
##                                              usesnplist[snp,3]-ydfile$pos[p]+1,
##                                              usesnplist[snp,3]-ydfile$pos[p]+1))
##       }
##       for(p in which(nkey==key)){
##         storealleles<-c(storealleles,substr(ndfile$seq[p],
##                                             usesnplist[snp,3]-ndfile$pos[p]+1,
##                                             usesnplist[snp,3]-ndfile$pos[p]+1))
##       }
##       
##       ortab<-outres[[1]][,,as.numeric(rownames(usesnplist)[snp])]
##       usebase<-names(sort(ortab,decreasing=T))[1:2]
##       
##       
##       allelecounts<-c(sum(storealleles==usebase[1]),sum(storealleles==usebase[2]))
##       if(allelecounts[2]>allelecounts[1]){allelecounts<-allelecounts[2:1]}
##       
##       HetSNPTableAABB[1,match(paste(allelecounts,collapse=":"),ACvec)]<-
##         HetSNPTableAABB[1,match(paste(allelecounts,collapse=":"),ACvec)]+1
##     }
##   }
## }
## 
## 
## 
## usesnplist<-snplistAB
## 
## #SNPnumbers[HSTrow]<-dim(usesnplist)[1]
## # now we are going to count up the numbers of duplicates that share the 
## # same allele and the numbers that do not.
## 
## # for each allele on the list
## for(snp in 1:dim(usesnplist)[1]){
##   if(100*trunc(snp/100)==snp){cat(snp,"\t")}
##   newwhich<-GRanges(usesnplist[snp,1], IRanges(usesnplist[snp,3], usesnplist[snp,3]))
##   
##   #these are the reads that were marked as duplicates 
##   ydfile<-scanBam(bamfile, 
##                   param=ScanBamParam(flag=scanBamFlag(isDuplicate=T),
##                                      simpleCigar=T,what=c("pos","mpos","seq"),
##                                      which=newwhich))[[1]]
##   
##   #these are the reads that were not
##   ndfile<-scanBam(bamfile,
##                   param=ScanBamParam(flag=scanBamFlag(isDuplicate=F),
##                                      simpleCigar=T,what=c("pos","mpos","seq"),
##                                      which=newwhich))[[1]]
##   
##   # Any fragement marked as a duplicate must be a duplicate of a fragment 
##   # that is not marked as a duplicate. We just want to keep a matched set 
##   # of duplicates and the fragments of which they are duplicates
##   
##   ykey<-paste(ydfile$pos,ydfile$mpos)
##   nkey<-paste(ndfile$pos,ndfile$mpos)
##   
##   usekey<-unique(ykey)
##   
##   # Now, assuming that we see some duplicates, we are going to go through them 
##   # and compare the alleles
##   if(length(usekey)>0){
##     for(key in usekey){
##       
##       storealleles<-NULL
##       for(p in which(ykey==key)){
##         storealleles<-c(storealleles, substr(ydfile$seq[p], 
##                                              usesnplist[snp,3]-ydfile$pos[p]+1,
##                                              usesnplist[snp,3]-ydfile$pos[p]+1))
##       }
##       for(p in which(nkey==key)){
##         storealleles<-c(storealleles,substr(ndfile$seq[p],
##                                             usesnplist[snp,3]-ndfile$pos[p]+1,
##                                             usesnplist[snp,3]-ndfile$pos[p]+1))
##       }
##       
##       ortab<-outres[[1]][,,as.numeric(rownames(usesnplist)[snp])]
##       usebase<-names(sort(ortab,decreasing=T))[1:2]
##       
##       
##       allelecounts<-c(sum(storealleles==usebase[1]),sum(storealleles==usebase[2]))
##       if(allelecounts[2]>allelecounts[1]){allelecounts<-allelecounts[2:1]}
##       
##       HetSNPTableAABB[2,match(paste(allelecounts,collapse=":"),ACvec)]<-
##         HetSNPTableAABB[2,match(paste(allelecounts,collapse=":"),ACvec)]+1
##     }
##   }
## }
## 
## HetSNPTableAABB<-cbind(c(10000,985),HetSNPTableAABB)
## 
## colnames(HetSNPTableAABB)<-c("NoSNPs", "AA", "AB", 
##                          "AAA", "AAB", 
##                          "AAAA", "AAAB", "AABB", 
##                          "AAAAA", "AAAAB", "AAABB", 
##                          "AAAAAA", "AAAAAB", "AAAABB", "AAABBB", 
##                          "AAAAAAA", "AAAAAAB", "AAAAABB", "AAAABBB")
## rownames(HetSNPTableAABB)<-c("AABB","AB")
## write.table(HetSNPTableAABB,file="HetSNPDupsCancer.txt",sep="\t")


###################################################
### code chunk number 51: LoadSNPsummaryTumour
###################################################
HetSNPTableAABB<-read.delim(file.path(EDpath,"Tumour","HetSNPDupsCancer.txt"),as.is=T)


###################################################
### code chunk number 52: ProcessTumourPatterns
###################################################
FDvecCancer<-rep(0,2)
for(j in 1:2){
  cat(j,"\t")
  vals<-rep(0,1000)
  myseq<-seq(0,0.2,length.out=1000)
  for(i in 1:1000){
    vals[i]<-propllik(HetSNPTableAABB[j,-1],myseq[i])
   }
  FDvecCancer[j]<-myseq[which.max(vals)]
}


###################################################
### code chunk number 53: TumourFradDup
###################################################
FDvecCancer


###################################################
### code chunk number 54: TumourConsequences
###################################################
# What were our PCR duplication rate estimates from Picard?
BasicAABB<-(TumourPicard[3,6]-TumourPicard[3,7])/(TumourPicard[3,3]-TumourPicard[3,7])
BasicAB<-(TumourPicard[2,6]-TumourPicard[2,7])/(TumourPicard[2,3]-TumourPicard[2,7])

# What are our corrected estimates?
CorAB<-(1-FDvecCancer[2])*BasicAB
CorAABB<-(1-FDvecCancer[1])*BasicAABB

# How much closer are they?

BasicAABB-BasicAB
CorAABB-CorAB
(BasicAABB-BasicAB)/(CorAABB-CorAB)


###################################################
### code chunk number 55: GetCEME
###################################################
CEME<-rep(22,0)
for(i in 1:22){
  R<-BasicEst[i]
  CEME[i]<-(readPairsExamined[i,1]--OpticalPairDuplicates[i,1])/(2*R)
}


###################################################
### code chunk number 56: GetCEbasic
###################################################
CEbasic<-rep(22,0)
for(i in 1:22){
  R<-BasicEst[i]
  N<-readPairsExamined[i,1]-OpticalPairDuplicates[i,1]
  startX<-CEME[i]
  CEbasic[i]<-optim(startX,libCompNewParam,R=R,N=N)$par
}


###################################################
### code chunk number 57: GetCEresid
###################################################
CEresid<-rep(22,0)
for(i in 1:22){
  R<-ResidEst[i]
  N<-readPairsExamined[i,1]-OpticalPairDuplicates[i,1]
  startX<-CEME[i]
  CEresid[i]<-optim(startX,libCompNewParam,R=R,N=N)$par
}


###################################################
### code chunk number 58: GetCEFDCor
###################################################
CEFDCor<-rep(22,0)
for(i in 1:22){
  R<-FDCorEst[i]
  N<-readPairsExamined[i,1]-OpticalPairDuplicates[i,1]
  startX<-CEME[i]
  CEFDCor[i]<-optim(startX,libCompNewParam,R=R,N=N)$par
}


###################################################
### code chunk number 59: GetCEFDCor2
###################################################
summary(CEFDCor/CEbasic)


###################################################
### code chunk number 60: Mito1
###################################################
lf<-read.delim(file.path(EDpath,"Picard",metfiles[3]),as.is=T)
lf<-lf[-(15:19),]


###################################################
### code chunk number 61: Mito2
###################################################
PicardMito<-cbind(sapply(lf[,1],substr,1,9),lf[,c(3,6,8)],
                  as.numeric(lf[,6])/as.numeric(lf[,3]),
                  as.numeric(lf[,5])/as.numeric(lf[,2]),ResDepth*(1-FDCorEst))


###################################################
### code chunk number 62: Mito3
###################################################
PicardMito<-cbind(PicardMito,as.numeric(PicardMito[,2])*(1-FDCorEst))


###################################################
### code chunk number 63: Mito4
###################################################
PicardMito<-cbind(PicardMito,as.numeric(PicardMito[,2])*400/(as.numeric(PicardMito[,7])*16569))


###################################################
### code chunk number 64: Mito5
###################################################
PicardMito<-cbind(PicardMito,
                  (as.numeric(PicardMito[,2])-as.numeric(PicardMito[,3]))
                  *400/(as.numeric(PicardMito[,7])*16569))
PicardMito<-cbind(PicardMito,as.numeric(PicardMito[,8])*400/(as.numeric(PicardMito[,7])*16569))
                  
colnames(PicardMito)<-c("Library","ReadPairs","Duplicates","PicardRate",
                        "BetterRate","UnpairedRate","ResCoverage",
                        "CorrectedDuplicateReadPairs","MTCNignoredups",
                        "MTCNremdups","MTCNcordups")


###################################################
### code chunk number 65: AddFig3
###################################################
par(mfrow=c(1,2))

plot(c(1,3),c(3,30),type="n",ylab="Library complexity estimate (billions)",xlab="",axes=F)
axis(2)
box()
for(i in seq(0,40,2)){
  rect(-1,i-1,5,i,border="grey95",col="grey95")
}
for(i in 1:22){
  tcol<-"black"
  if(WeaverSuppTable1[i,12]=="Blood"){tcol<-"red"}
  lines(1:3,c(CEbasic[i],CEresid[i],CEFDCor[i])/(10^9),type="b",pch=16,col=tcol)
}
axis(1,at=1:3,labels=c("Basic\nduplicate\nestimate",
                       "Residual\nduplicate\nestimate",
                       "Corrected\nestimate"),las=3)

plot(c(1,3),c(60,800),type="n",ylab="MtDNA copy number",xlab="",axes=F)
axis(2)
box()
for(i in seq(0,800,200)){
  rect(-1,i-100,5,i,border="grey95",col="grey95")
}
for(i in 1:22){
  tcol<-"black"
  if(WeaverSuppTable1[i,12]=="Blood"){tcol<-"red"}
  lines(1:3,(as.numeric(PicardMito[i,c(9,11,10)])),type="b",pch=16,col=tcol)
}
axis(1,at=1:3,labels=c("No\nduplicates\nremoved","New\nestimate","All\nduplicates\nremoved"),las=3)


###################################################
### code chunk number 66: Session Info
###################################################
toLatex(sessionInfo())


