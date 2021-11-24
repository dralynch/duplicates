genprobsmat <-
  function(N,FD,PL,CL,reso=1001){
    
    PD<-1-FD
    
    index<-1
    usePL<-PL[[N]]
    
    params<-matrix(NA,nrow=dim(usePL)[1],ncol=3)
    
    
    params[,2]<-apply(usePL,1,function(x){sum(x>0)-1})
    params[,1]<-N-params[,2]-1
    params[,3]<-CL[[N]]
    
    partitionprobs<-matrix(NA,ncol=dim(usePL)[1],nrow=length(FD))
    
    ## Now calculate the probability of each partition
    for(i in 1:dim(usePL)[1]){
      partitionprobs[,i]<-params[i,3]*(PD^params[i,1])*(FD^params[i,2])
    }
    
    
    ## noB = number of B alleles
    noB<-0:floor(N/2)
    ## probnoB = probability of seeing that number of B alleles
    probnoB<-matrix(0,ncol=length(noB),nrow=length(FD))
    
    nop<-dim(usePL)[2]
    binky<-as.matrix(sapply(0:((2^(nop-1))-1),function(x){as.integer(intToBits(x))})[1:nop,])
    
    for(part in 1:dim(usePL)[1]){
      ## nop = number of partitions
      mytab<-(rev(usePL[part,])%*%binky)
      mytab[mytab>nop/2]<-nop-mytab[mytab>nop/2]
      mytab<-table(mytab)
      for(k in 1:length(mytab)){
        probnoB[,as.integer(names(mytab)[k])+1]<-probnoB[,as.integer(names(mytab)[k])+1]+mytab[k]*partitionprobs[,part]/(2^(nop-1))        
      }
      
    }
    
    return(probnoB)
  }


processduptable<-function(HST,PL=NULL,CL=NULL,reso=1001){
  
  # Check whether there is a partition list
  if(is.null(PL)){
    PL<-genPL(nchar(colnames(HST)[dim(HST)[2]-1]))
    CL<-genCoefs(PL)
  }else{
    # Check whether the existing partition list will suffice
    if(nchar(colnames(HST)[dim(HST)[2]-1])>length(PL)){
      message("NEED TO GENERATE A LARGER PARTITION LIST")
      PL<-genPL(nchar(colnames(HST)[dim(HST)[2]-1]))
      CL<-genCoefs(PL)
    }
    if(is.null(CL)){
      CL<-genCoefs(PL)
    }
  }
  
  FDvec<-rep(0,dim(HST)[1])
  
  myseq<-seq(0,1,length.out=reso)[-1]
  probmat<-matrix(0,ncol=(dim(HST)[2]-4),nrow=(reso-1))
  index<-0
  for(i in 2:nchar(colnames(HST)[dim(HST)[2]-1])){
    #message(i)
    tempmat<-genprobsmat(i,myseq,PL,CL,reso=reso)
    probmat[,index+(1:(dim(tempmat)[2]))]<-tempmat
    index<-index+dim(tempmat)[2]
  }
  
  myvecs<-as.matrix(HST[,-c(1:3,dim(HST)[2]),drop=FALSE]) %*% t(log(probmat))
  FDvec<-myseq[apply(myvecs,1,which.max)]
    
  outtab<-cbind(100*HST[,3]/HST[,2],FDvec,100*HST[,3]/HST[,2]*(1-FDvec),100*HST[,3]/HST[,2]*(FDvec))
  colnames(outtab)<-c("ObservedDupRate","PropFragDups","AdjustedDupRate","FragDupRate")
  rownames(outtab)<-rownames(HST)
  return(outtab)
}
