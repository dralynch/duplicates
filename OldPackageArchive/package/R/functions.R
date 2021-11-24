genCoefs <-
  function(PL){
    
    CoefList<-list(1)
    for(i in 2:length(PL)){
      oldPL<-PL[[i-1]]
      newPL<-PL[[i]]
      oldlist<-CoefList[[i-1]]
      newlist<-rep(0,dim(newPL)[1])
      
      for(k in 1:dim(oldPL)[1]){
        #cat("k=",k,"\n")
        oldcoef<-oldlist[k]
        oldrow<-oldPL[k,]
        lor<-sum(oldrow>0)
        for(l in 1:lor){
          #cat("l=",l,"\n")        
          newrow<-oldrow
          newrow[l]<-newrow[l]+1
          newrow<-c(newrow,0)
          newrow<-sort(newrow,decreasing=T)
          for(m in 1:dim(newPL)[1]){
            if(all(newPL[m,]==newrow)){newlist[m]<-newlist[m]+oldcoef/lor}
          }
        }
        newrow<-sort(c(oldrow,1),decreasing=T)
        for(m in 1:dim(newPL)[1]){
          if(all(newPL[m,]==newrow)){newlist[m]<-newlist[m]+oldcoef}
        }
      }
      CoefList[[i]]<-newlist
    }
    return(CoefList)
  }

genPL <- 
  function(LPL=20){
    if(LPL<4){LPL<-4}
    ## Hard-code the first two entries
    partslist<-list(matrix(1,ncol=1,nrow=1),matrix(c(2,0,1,1),ncol=2,nrow=2,byrow=T))
    index<-3
    for(i in 3:LPL){
      Lcol<-(i-1):1
      newel<-matrix(c(i,rep(0,i-1)),nrow=1)
      for(j in Lcol){
        temp<-as.matrix(partslist[[i-j]])
        newel<-rbind(newel,cbind(j,temp,matrix(0,ncol=(i-(dim(temp)[2])-1),nrow=dim(temp)[1])))
      }
      newel<-newel[newel[,1]>=apply(newel,1,max),]
      colnames(newel)<-1:(dim(newel)[2])
      partslist[[index]]<-newel
      index<-index+1
    }
    return(partslist)
  }


genprobs <-
  function(N,FD,fullout=F,quickcoef=F){
    ## function to work out the AAAA,AAAB,AABB etc. probabilities given
    ## the Number of duplicates and the fragmentation duplication rate
    
    ## For neatness, define the PCR duplication rate
    PD<-1-FD
    
    ## identify the list of possible partitions into distinct molecules
    ## given the number of duplicates
    
    partitions<-list(NULL)
    index<-1
    usepl<-partitionlist[[N]]
    for(i in 1:(dim(usepl)[1])){
      temp<-as.vector(usepl[i,])
      temp<-temp[temp!=0]
      partitions[[i]]<-temp
      index<-index+1
    }
    
    ## for each partition work out the number of PCR duplicates,
    ## the number of fragmentation duplicates, and the
    ## contribution to the binomial coefficient
    ## This latter term can be read from a pre-compiled table, or 
    ## generated on the fly if preferred 
    params<-matrix(NA,nrow=length(partitions),ncol=3)
    
    for(i in 1:length(partitions)){
      params[i,2]<-length(partitions[[i]])-1
      params[i,1]<-N-params[i,2]-1
      if(quickcoef){
        params[i,3]<-factorial(length(partitions[[i]]))/prod(factorial(table(partitions[[i]])))
      }
    }
    if(!quickcoef){
      params[,3]<-CoefList[[N]]
    }
    ## Now calculate the probability of each partition
    partitionprobs<-rep(NA,length(partitions))
    for(i in 1:length(partitions)){
      partitionprobs[i]<-params[i,3]*(PD^params[i,1])*(FD^params[i,2])
    }
    
    
    ## Next, we take each partition of duplicates into distinct molecules
    ## and consider the possible observed alleles that could be seen
    ## assuming that the patient is diploid and heterozygous at the site
    
    ## For the purposes of tallying, we always define the A allele as being the 
    ## (equal-)most frequent allele. Therefore the B allele count can be at most
    ## half of the number of duplicates
    
    ## noB = number of B alleles
    noB<-0:floor(N/2)
    ## probnoB = probability of seeing that number of B alleles
    probnoB<-rep(0,length(noB))
    
    for(part in 1:length(partitions)){
      ## nop = number of partitions
      nop<-length(partitions[[part]])
      ## binky is the matrix where each column represents an equally likely
      ## allocation of alleles to molecules
      ## binky = 'binary key'
      ## this is not a harry potter reference
      ## it is a discworld reference 
      ## that is much more reasonable
      binky<-as.matrix(sapply(0:((2^(nop-1))-1),intToBits)[1:nop,])
      for(ABsplit in 1:2^(nop-1)){
        ABalloc<-binky[,ABsplit]
        ABalloc<-ABalloc[length(ABalloc):1]
        ## tempnoB is the temporary B count for this allocation of alleles
        tempnoB<-sum(partitions[[part]]*as.integer(ABalloc))
        if(tempnoB>floor(N/2)){tempnoB<-N-tempnoB}
        probnoB[tempnoB+1]<-probnoB[tempnoB+1]+partitionprobs[part]/(2^(nop-1))
      }  
    }
    
    retval<-probnoB
    if(fullout){retval<-list(partitions,params,partitionprobs,probnoB)}
    return(retval)
  }


libCompNewParam <-
  function(X,R,N){
    (log(N*(R-1)+X)-log(X)+N/X)^2
  }

propllik <-
  function(x,FD){
    ## function that takes the empirical observations of allele patterns,
    ## and calculates the log-likelihood of the parameter FD, the 
    ## proportion of duplicates arising from fragmentation of distinct 
    ## molecules, via calls to genprobs.
    mysum<-0
    index<-0
    for(i in 1:6){    
      myp<-genprobs(i+1,FD,fullout=F)
      mysum<-mysum+sum(x[(index+1):(index+length(myp))]*log(myp))
      index<-index+length(myp)
    }
    return(mysum)
  }

secondbiggest <-
  function(x){
    sort(x, decreasing = T)[2]
  }
