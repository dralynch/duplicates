
processBAMfordups<-function(bamfilelist,path,snplist,chromcol=1,poscol=3,maxdups=7,AFcutoff=0.4,maxdepth=200,mindepth=0){
  snplist<-snplist[which(snplist[,chromcol]!="chrY"),]
  
  
  #produce the list of patterns being considered
  ACvec<-NULL
  HSTCN<-c("NoSNPs","NoReads","NoDups")
  for(i in 2:maxdups){
    for(j in 0:floor(i/2)){
      ACvec<-c(ACvec,paste(i-j,j,sep=":",collapse=""))
      HSTCN<-c(HSTCN,paste(c(rep("A",i-j),rep("B",j)),sep="",collapse=""))
    }
  }
  HSTCN<-c(HSTCN,"NoN")
  
  HetSNPTable<-matrix(0,nrow=length(bamfilelist),ncol=length(ACvec))
  SNPnumbers<-rep(0,length(bamfilelist))
  
  Readnumbers<-rep(0,length(bamfilelist))
  Dupnumbers<-rep(0,length(bamfilelist))
  Nnumbers<-rep(0,length(bamfilelist))
  
  #cycle through bam files
  for(HSTrow in 1:length(bamfilelist)){
    sample<-bamfilelist[HSTrow]
    message(sample)
    bamfile=paste(path,sample,sep="")
    # look to see which of the candidate SNPs is heterozygous in this sample.
    # Note that this is based on all the reads not just the duplicate fragments,
    # so should not greatly bias matters.
    fls <- PileupFiles(bamfile) 
    which<-GRanges(snplist[,chromcol],IRanges(snplist[,poscol],snplist[,poscol]))
    PUP <- ApplyPileupsParam(which=which, yieldSize=1000000L, yieldBy="position", what="seq",maxDepth=maxdepth,minDepth=mindepth)
    outres <- applyPileups(fls,(function(x){x[["seq"]]}),param=PUP)
    snpids<-which(apply((outres)[[1]],3,secondbiggest)/apply((outres)[[1]],3,sum)>AFcutoff)
    usesnplist<-snplist[ snpids,]
    
    SNPnumbers[HSTrow]<-dim(usesnplist)[1]
    
    # now we are going to count up the numbers of duplicates that share the
    # same allele and the numbers that do not.
    
    
    # for each allele on the list
    for(snp in 1:dim(usesnplist)[1]){
      
      ortab<-outres[[1]][,,snpids[snp]]
      usebase<-names(sort(ortab,decreasing=T))[1:2]          
      
      #cat(snp,"\t")
      newwhich<-GRanges(usesnplist[snp,chromcol], IRanges(usesnplist[snp,poscol], usesnplist[snp,poscol]))
      
      zdfile<-scanBam(bamfile, param=ScanBamParam(flag=scanBamFlag(isProperPair=T),simpleCigar=T,what=c("pos","seq","flag","isize"),
                                                  which=newwhich))[[1]]
      
      
      
      
      
      
      
      
      
      
      
      
      zdkey<-paste(zdfile$pos,zdfile$isize)
      
      
      if(length(zdfile$flag)>0){        
        ykey<-zdkey[which(sapply(zdfile$flag,function(x){as.logical(intToBits(x)[11])}))]
        ykey<-ykey[ykey %in% zdkey[-(match(ykey,zdkey))]]
        
        Readnumbers[HSTrow]<-Readnumbers[HSTrow]+length(zdkey)
        Dupnumbers[HSTrow]<-Dupnumbers[HSTrow]+length(ykey)
        
        ykey<-unique(ykey)
        
        # Now, assuming that we see some duplicates, we are going to go through them 
        # and compare the alleles
        if(length(ykey)>0){
          for(key in ykey){
            
            storealleles<-NULL
            for(p in which(zdkey==key)){
              storealleles<-c(storealleles, substr(zdfile$seq[p], 
                                                   usesnplist[snp,poscol]-zdfile$pos[p]+1,
                                                   usesnplist[snp,poscol]-zdfile$pos[p]+1))
            }
            Nnumbers[HSTrow]<-Nnumbers[HSTrow]+sum(!(storealleles %in% usebase))
            
            allelecounts<-c(sum(storealleles==usebase[1]),sum(storealleles==usebase[2]))
            if(allelecounts[2]>allelecounts[1]){allelecounts<-allelecounts[2:1]}
            
            HetSNPTable[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]<-HetSNPTable[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]+1
            if(is.na(match(paste(allelecounts,collapse=":"),ACvec))){
              Readnumbers[HSTrow]<-Readnumbers[HSTrow]-sum(zdkey==key)
              Dupnumbers[HSTrow]<-Dupnumbers[HSTrow]-sum(zdkey==key)+1
              Nnumbers[HSTrow]<-Nnumbers[HSTrow]-sum(!(storealleles %in% usebase))
            }
          }# duplicates key cycle
        }#duplicates key exists
      }
    }# snp cycle
    gc()
  }# bamfile cycle
  
  
  HetSNPTable<-cbind(SNPnumbers,Readnumbers,Dupnumbers,HetSNPTable,Nnumbers)
  colnames(HetSNPTable)<-HSTCN
  rownames(HetSNPTable)<-bamfilelist
  return(HetSNPTable)
}
processBAMbylane<-function(bamfilelist,path,snplist,lanelist,chromcol=1,poscol=3,maxdups=7,AFcutoff=0.4,maxdepth=200,mindepth=0){
  snplist<-snplist[which(snplist[,chromcol]!="chrY"),]
  if(length(bamfilelist)>1){stop("bamfilelist must be of length 1")}
  
  #produce the list of patterns being considered
  ACvec<-NULL
  HSTCN<-c("NoSNPs","NoReads","NoDups")
  for(i in 2:maxdups){
    for(j in 0:floor(i/2)){
      ACvec<-c(ACvec,paste(i-j,j,sep=":",collapse=""))
      HSTCN<-c(HSTCN,paste(c(rep("A",i-j),rep("B",j)),sep="",collapse=""))
    }
  }
  HSTCN<-c(HSTCN,"NoN")
  
  HetSNPTable<-matrix(0,nrow=(length(lanelist)),ncol=length(ACvec))
  
  
  Readnumbers<-rep(0,length(lanelist))
  Dupnumbers<-rep(0,length(lanelist))
  Nnumbers<-rep(0,length(bamfilelist))
  
  
  
  sample<-bamfilelist
  message(sample)
  bamfile=paste(path,sample,sep="")
  # look to see which of the candidate SNPs is heterozygous in this sample.
  # Note that this is based on all the reads not just the duplicate fragments,
  # so should not greatly bias matters.
  fls <- PileupFiles(bamfile) 
  which<-GRanges(snplist[,chromcol],IRanges(snplist[,poscol],snplist[,poscol]))
  PUP <- ApplyPileupsParam(which=which, yieldSize=1000000L, yieldBy="position", what="seq",maxDepth=maxdepth,minDepth=mindepth)
  outres <- applyPileups(fls,(function(x){x[["seq"]]}),param=PUP)
  snpids<-which(apply((outres)[[1]],3,secondbiggest)/apply((outres)[[1]],3,sum)>AFcutoff)
  usesnplist<-snplist[ snpids,]
  
  SNPnumbers<-rep(dim(usesnplist)[1],length(lanelist))
  
  # now we are going to count up the numbers of duplicates that share the
  # same allele and the numbers that do not.
  
  
  # for each allele on the list
  for(snp in 1:dim(usesnplist)[1]){
    cat(snp,"\t")
    ortab<-outres[[1]][,,snpids[snp]]
    usebase<-names(sort(ortab,decreasing=T))[1:2]          
    
    
    newwhich<-GRanges(usesnplist[snp,chromcol], IRanges(usesnplist[snp,poscol], usesnplist[snp,poscol]))
    
    zdfile<-scanBam(bamfile, param=ScanBamParam(flag=scanBamFlag(isProperPair=T),simpleCigar=T,what=c("qname","pos","seq","flag","isize"),
                                                which=newwhich))[[1]]
    
    for(mylane in 1:length(lanelist)){
      
      zdfile1<-zdfile
      
      zdfile1$seq<-zdfile1$seq[sapply(zdfile$qname,substr,1,13) %in% lanelist[[mylane]]]
      zdfile1$pos<-zdfile1$pos[sapply(zdfile$qname,substr,1,13) %in% lanelist[[mylane]]]
      #zdfile1$mpos<-zdfile1$mpos[sapply(zdfile$qname,substr,1,13) %in% lanelist[[mylane]]]
      zdfile1$flag<-zdfile1$flag[sapply(zdfile$qname,substr,1,13) %in% lanelist[[mylane]]]
      zdfile1$isize<-zdfile1$isize[sapply(zdfile$qname,substr,1,13) %in% lanelist[[mylane]]]
      zdfile1$qname<-zdfile1$qname[sapply(zdfile$qname,substr,1,13) %in% lanelist[[mylane]]]
      
      zdkey<-paste(zdfile1$pos,zdfile1$isize)
      
      
      if(length(zdfile1$flag)>0){
        ykey<-zdkey[which(sapply(zdfile1$flag,function(x){as.logical(intToBits(x)[11])}))]
        ykey<-ykey[ykey %in% zdkey[-(match(ykey,zdkey))]]
        
        Readnumbers[mylane]<-Readnumbers[mylane]+length(zdkey)
        Dupnumbers[mylane]<-Dupnumbers[mylane]+length(ykey)
        
        ykey<-unique(ykey)
        
        # Now, assuming that we see some duplicates, we are going to go through them 
        # and compare the alleles
        if(length(ykey)>0){
          for(key in ykey){
            
            storealleles<-NULL
            for(p in which(zdkey==key)){
              storealleles<-c(storealleles, substr(zdfile1$seq[p], 
                                                   usesnplist[snp,poscol]-zdfile1$pos[p]+1,
                                                   usesnplist[snp,poscol]-zdfile1$pos[p]+1))
            }
            Nnumbers[mylane]<-Nnumbers[mylane]+sum(!(storealleles %in% usebase))
            
            allelecounts<-c(sum(storealleles==usebase[1]),sum(storealleles==usebase[2]))
            if(allelecounts[2]>allelecounts[1]){allelecounts<-allelecounts[2:1]}
            
            HetSNPTable[mylane,match(paste(allelecounts,collapse=":"),ACvec)]<-HetSNPTable[mylane,match(paste(allelecounts,collapse=":"),ACvec)]+1  
            if(is.na(match(paste(allelecounts,collapse=":"),ACvec))){
              Readnumbers[mylane]<-Readnumbers[mylane]-sum(zdkey==key)
              Dupnumbers[mylane]<-Dupnumbers[mylane]-sum(zdkey==key)+1
              Nnumbers[mylane]<-Nnumbers[mylane]-sum(!(storealleles %in% usebase))
            }      
          }# duplicates key cycle
        }#duplicates key exists
      }
    }# lane cycle
  }# snp cycle
  
  
  HetSNPTable<-cbind(SNPnumbers,Readnumbers,Dupnumbers,HetSNPTable,Nnumbers)
  colnames(HetSNPTable)<-HSTCN
  rownames(HetSNPTable)<-1:length(lanelist)
  return(HetSNPTable)
}

processBAMSE<-function(bamfilelist,path,snplist,chromcol=1,poscol=3,maxdups=7,AFcutoff=0.4,maxdepth=200,mindepth=0){
  snplist<-snplist[which(snplist[,chromcol]!="chrY"),]
  
  
  #produce the list of patterns being considered
  ACvec<-NULL
  HSTCN<-c("NoSNPs","NoReads","NoDups")
  for(i in 2:maxdups){
    for(j in 0:floor(i/2)){
      ACvec<-c(ACvec,paste(i-j,j,sep=":",collapse=""))
      HSTCN<-c(HSTCN,paste(c(rep("A",i-j),rep("B",j)),sep="",collapse=""))
    }
  }
  HSTCN<-c(HSTCN,"NoN")
  
  HetSNPTable<-matrix(0,nrow=length(bamfilelist),ncol=length(ACvec))
  SNPnumbers<-rep(0,length(bamfilelist))
  
  Readnumbers<-rep(0,length(bamfilelist))
  Dupnumbers<-rep(0,length(bamfilelist))
  Nnumbers<-rep(0,length(bamfilelist))
  
  #cycle through bam files
  for(HSTrow in 1:length(bamfilelist)){
    sample<-bamfilelist[HSTrow]
    message(sample)
    bamfile=paste(path,sample,sep="")
    # look to see which of the candidate SNPs is heterozygous in this sample.
    # Note that this is based on all the reads not just the duplicate fragments,
    # so should not greatly bias matters.
    fls <- PileupFiles(bamfile) 
    which<-GRanges(snplist[,chromcol],IRanges(snplist[,poscol],snplist[,poscol]))
    PUP <- ApplyPileupsParam(which=which, yieldSize=1000000L, yieldBy="position", what="seq",maxDepth=maxdepth,minDepth=mindepth)
    outres <- applyPileups(fls,(function(x){x[["seq"]]}),param=PUP)
    snpids<-which(apply((outres)[[1]],3,secondbiggest)/apply((outres)[[1]],3,sum)>AFcutoff)
    usesnplist<-snplist[ snpids,]
    
    SNPnumbers[HSTrow]<-dim(usesnplist)[1]
    
    # now we are going to count up the numbers of duplicates that share the
    # same allele and the numbers that do not.
    
    
    # for each allele on the list
    for(snp in 1:dim(usesnplist)[1]){
      
      ortab<-outres[[1]][,,snpids[snp]]
      usebase<-names(sort(ortab,decreasing=T))[1:2]          
      
      #cat(snp,"\t")
      newwhich<-GRanges(usesnplist[snp,chromcol], IRanges(usesnplist[snp,poscol], usesnplist[snp,poscol]))
      
      zdfile<-scanBam(bamfile, param=ScanBamParam(flag=scanBamFlag(isProperPair=T),simpleCigar=T,what=c("pos","seq","flag"),
                                                  which=newwhich))[[1]]
      
      zdkey<-paste(zdfile$pos)
      
      
      if(length(zdfile$flag)>0){        
        #ykey<-zdkey[which(sapply(zdfile$flag,function(x){as.logical(intToBits(x)[11])}))]
        #ykey<-ykey[ykey %in% zdkey[-(match(ykey,zdkey))]]
        
        ykey<-zdkey[duplicated(zdkey)]
        
        Readnumbers[HSTrow]<-Readnumbers[HSTrow]+length(zdkey)
        Dupnumbers[HSTrow]<-Dupnumbers[HSTrow]+length(ykey)
        
        ykey<-unique(ykey)
        
        # Now, assuming that we see some duplicates, we are going to go through them 
        # and compare the alleles
        if(length(ykey)>0){
          for(key in ykey){
            
            storealleles<-NULL
            for(p in which(zdkey==key)){
              storealleles<-c(storealleles, substr(zdfile$seq[p], 
                                                   usesnplist[snp,poscol]-zdfile$pos[p]+1,
                                                   usesnplist[snp,poscol]-zdfile$pos[p]+1))
            }
            Nnumbers[HSTrow]<-Nnumbers[HSTrow]+sum(!(storealleles %in% usebase))
            
            allelecounts<-c(sum(storealleles==usebase[1]),sum(storealleles==usebase[2]))
            if(allelecounts[2]>allelecounts[1]){allelecounts<-allelecounts[2:1]}
            
            HetSNPTable[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]<-HetSNPTable[HSTrow,match(paste(allelecounts,collapse=":"),ACvec)]+1
            if(is.na(match(paste(allelecounts,collapse=":"),ACvec))){
              Readnumbers[HSTrow]<-Readnumbers[HSTrow]-sum(zdkey==key)
              Dupnumbers[HSTrow]<-Dupnumbers[HSTrow]-sum(zdkey==key)+1
              Nnumbers[HSTrow]<-Nnumbers[HSTrow]-sum(!(storealleles %in% usebase))
            }
          }# duplicates key cycle
        }#duplicates key exists
      }
    }# snp cycle
    gc()
  }# bamfile cycle
  HetSNPTable<-cbind(SNPnumbers,Readnumbers,Dupnumbers,HetSNPTable,Nnumbers)
  colnames(HetSNPTable)<-HSTCN
  rownames(HetSNPTable)<-bamfilelist
  return(HetSNPTable)
  
}
