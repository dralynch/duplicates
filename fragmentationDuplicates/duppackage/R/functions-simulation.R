
simreads<-function(readlength=100,depth=50,pcrdup=0.05,ISmean=310,ISsd=7,realISD=NULL,genlength=999999,spacing=1000,lowerlength=200,upperlength=400,maxM=NULL){
  
  numfrags<-depth*(1-pcrdup)*genlength/(2*readlength)
  
  hetsnplocs<-seq(spacing,genlength,spacing)
  
  fragstarts<-sample(genlength,numfrags,replace=T)
  if(is.null(realISD)){fraglengths<-round(rnorm(numfrags,ISmean,ISsd^2))}
  else{fraglengths<-sample(realISD[,1],numfrags,prob=realISD[,2],replace=T)}
  
  
  
  keeplength<-which((fraglengths>=lowerlength)&(fraglengths<=upperlength))
  
  fragstarts<-fragstarts[keeplength]
  fraglengths<-fraglengths[keeplength]
  
  fragends<-fraglengths+fragstarts
  
  keepsnp<-which(((spacing - (fragstarts %% spacing))<readlength)|(((fragends %% spacing))<readlength))
  
  fragstarts<-fragstarts[keepsnp]
  fraglengths<-fraglengths[keepsnp]
  
  snpstate<-sample(c(0,1),length(fragstarts),replace=T)
  
  pcrdupnos<-rpois(length(fragstarts),pcrdup/(1-pcrdup))
  
  simdata<-cbind(fragstarts,fraglengths,snpstate)
  
  simdata<-simdata[rep(1:length(fragstarts),pcrdupnos+1),]
  #for(i in which(pcrdupnos>0)){
  #  for(j in 1:pcrdupnos[i]){
  #    simdata<-rbind(simdata,simdata[i,])
  #  }
  #}
  
  simkey<-paste(simdata[,1],simdata[,2])
  
  maxdups<-max(table(simkey))
  if(!is.null(maxM)){
    if(maxM>maxdups){
      maxdups<-maxM
    }
  }
  
  ACvec<-NULL
  HSTCN<-c("NoSNPs","NoReads","NoDups")
  for(i in 2:maxdups){
    for(j in 0:floor(i/2)){
      ACvec<-c(ACvec,paste(i-j,j,sep=":",collapse=""))
      HSTCN<-c(HSTCN,paste(c(rep("A",i-j),rep("B",j)),sep="",collapse=""))
    }
  }
  HSTCN<-c(HSTCN,"NoN")
  
  SimSNPTable<-matrix(0,nrow=1,ncol=length(ACvec))
  SNPnumbers<-length(hetsnplocs)
  
  Readnumbers<-length(simkey)
  Dupnumbers<-sum(duplicated(simkey))
  Nnumbers<-0
  
  dupkeylist<-unique(simkey[duplicated(simkey)])
  
  for(key in dupkeylist){
    tmpstates<-simdata[which(simkey==key),3]
    tmpstr<-paste(sort(c(sum(tmpstates==0),sum(tmpstates==1)),decreasing=T),collapse=":")
    SimSNPTable[1,match(tmpstr,ACvec)]<-SimSNPTable[1,match(tmpstr,ACvec)]+1
  }
  SimSNPTable<-cbind(SNPnumbers,Readnumbers,Dupnumbers,SimSNPTable,Nnumbers)
  colnames(SimSNPTable)<-HSTCN
  
  return(SimSNPTable)
}


simMTreads<-function(readlength=100,depthpercopy=25,numcopies,pcrdup=0.05,ISmean=310,ISsd=7,realISD=NULL,genlength=16500,spacing=1000,lowerlength=100,upperlength=500){
  
  #numfrags<-depthpercopy*numcopies*(1-pcrdup)*genlength/(2*readlength)
  numfrags<-depthpercopy*numcopies*genlength/(2*readlength)
  
  fragstarts<-sample(genlength,numfrags,replace=T)
  if(is.null(realISD)){fraglengths<-round(rnorm(numfrags,ISmean,ISsd^2))}
  else{fraglengths<-sample(realISD[,1],numfrags,prob=realISD[,2],replace=T)}
  
  keeplength<-which((fraglengths>=lowerlength)&(fraglengths<=upperlength))
  
  fragstarts<-fragstarts[keeplength]
  fraglengths<-fraglengths[keeplength]
  
  
  fragends<-fraglengths+fragstarts
  
  pcrdupnos<-rpois(length(fragstarts),pcrdup/(1-pcrdup))
  
  simdata<-cbind(fragstarts,fraglengths)
  
  simdata<-simdata[rep(1:length(fragstarts),pcrdupnos+1),]
  
  #for(i in which(pcrdupnos>0)){
  #  for(j in 1:pcrdupnos[i]){
  #   simdata<-rbind(simdata,simdata[i,])
  #  }
  #}
  
  simkey<-paste(simdata[,1],simdata[,2])
  
  
  Readnumbers<-length(simkey)
  Dupnumbers<-sum(duplicated(simkey))
  
  
  resout<-c(Readnumbers,Readnumbers*(1-pcrdup),Readnumbers-Dupnumbers,Readnumbers*2*readlength/genlength,Readnumbers*2*readlength*(1-pcrdup)/genlength,(Readnumbers-Dupnumbers)*2*readlength/genlength,Readnumbers*2*readlength/(genlength*depthpercopy),Readnumbers*2*readlength*(1-pcrdup)/(genlength*depthpercopy),(Readnumbers-Dupnumbers)*2*readlength/(genlength*depthpercopy))
  
  names(resout)<-c("ReadNos","GenuineReadNos","FilteredReadNos","ReadDepth","GenuineReadDepth","FilteredReadDepth","CopyNo","GenuineCopyNo","FilteredCopyNo")
  
  return(format(resout,scientific=F))
}
