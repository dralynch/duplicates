
partitionToCoef<-function(x){
  x<-x[x!=0]
  return(factorial(length(x))/prod(factorial(table(x))))
}

genCoefs <-
  function(PL){
    PTM<-proc.time()
    CoefList<-list(1)
    for(i in 2:length(PL)){
      message(i)
      CoefList[[i]]<-apply(PL[[i]],1,partitionToCoef)
      message(paste(proc.time() - PTM,collapse=" : "))
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

secondbiggest <-
  function(y){
    sort(y, decreasing = T)[2]
  }
