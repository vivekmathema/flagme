

betweenAlignment<-function(pD,cAList,pAList,impList,filterMin=3,gap=0.7,D=10,usePeaks=TRUE,df=30,verbose=TRUE) {
  n<-length(pAList)
  if(length(filterMin)==1) filterMin<-rep(filterMin,n)
  pkd<-vector("list",n)
  pkr<-vector("list",n)
  filtind<-vector("list",n)  # list of filtered indices
  newind<-vector("list",n)
  for(g in 1:n) {
    pa<-pAList[[g]]
	nm<-length(pa@merges)
    runs<-pa@merges[[nm]]$runs
    ind<-pa@merges[[nm]]$ind
	keep<-rowSums(!is.na(ind))>filterMin[g]
	if (verbose)
	  cat("[betweenAlignment]", names(pAList)[g], ": Removing", nrow(ind)-sum(keep),"peaks.\n")
    ind<-ind[keep,]
	newind[[g]]<-lapply(impList[[g]],FUN=function(u,k) u[k,],k=keep)
	filtind[[g]]<-ind
    rt<-matrix(NA,nrow(ind),ncol(ind))
    mz<-pD@mz
    peaksdata<-matrix(0,nrow=length(pD@mz),ncol=nrow(ind))
    for(j in 1:ncol(ind)) {
	  if (usePeaks) {
	    cur.ds<-pD@peaksdata[[runs[j]]]
		cur.rt<-pD@peaksrt
	  } else {
	    cur.ds<-pD@rawdata[[runs[j]]]
		cur.rt<-pD@rawrt
	  }
	  for(i in 1:nrow(ind)) {
		if(!is.na(ind[i,j])) {
		  peaksdata[,i]<-peaksdata[,i]+cur.ds[,ind[i,j]]
		  rt[i,j]<-cur.rt[[runs[j]]][ind[i,j]]
		}
	  }
    }
	pkd[[g]]<-peaksdata
	if (!usePeaks)
	  pkd[[g]]<-pkd[[g]] / outer(rep(1,nrow(peaksdata)),rowSums(!is.na(ind))) # average over the number of samples
	pkr[[g]]<-apply(rt,1,median,na.rm=TRUE)
  }
  #wRA<-new("peaksDataset",mz=mz,rawdata=NULL,files=names(cAList),rawrt=NULL,peaksdata=pkd,peaksrt=pkr,filtind=filtind)
  wRA<-new("peaksDataset",mz=mz,rawdata=list(NULL),files=names(cAList),rawrt=list(NULL),peaksdata=pkd,peaksrt=pkr)
  #return(wRA)
  #filtind<-
  cA<-clusterAlignment(wRA,1:n,gap=gap,D=D,df=df,verbose=verbose)
  pA<-progressiveAlignment(wRA,cA,gap=gap,D=D,df=df,verbose=verbose)
  nm<-length(pA@merges)
  ind<-pA@merges[[nm]]$ind
  groups<-pA@merges[[nm]]$runs
  full.runs<-NULL
  full.groups<-NULL
  full.ind<-matrix(NA,nrow(ind),length(pD@rawdata))
  full.newind<-vector("list",3)
  for(i in 1:length(full.newind))
    full.newind[[i]]<-matrix(NA,nrow(ind),length(pD@rawdata))
  names(full.newind)<-c("apex","start","end")
  
  col<-1
  for(i in 1:length(groups)) {
    g<-groups[i]
    n<-length(pAList[[g]]@merges)
    full.runs<-c(full.runs,pAList[[g]]@merges[[n]]$runs)
	nr<-length(pAList[[g]]@merges[[n]]$runs)
    full.groups<-c(full.groups,rep(wRA@files[g],nr))
	this.ind<-filtind[[g]]
	rw<-which(!is.na(ind[,i]))
	cl<-col:(col+ncol(this.ind)-1)
	full.ind[rw,cl]<-this.ind
    if(!is.null(newind[[g]]$apex)) {
	  for(j in 1:length(full.newind))
	    full.newind[[j]]<-newind[[g]][[j]]
	}
	col<-col+ncol(this.ind)
  }
  new("betweenAlignment",mergedPeaksDataset=wRA,ind=full.ind,imputeind=full.newind,runs=full.runs,groups=full.groups,cA=cA,pA=pA,filtind=filtind,newind=newind)
}

setMethod("show","betweenAlignment",
function(object) {
   cat("An object of class \"",class(object),"\"\n",sep="")
   cat(length(object@mergedPeaksDataset@peaksrt),"groups:",sapply(object@mergedPeaksDataset@peaksrt,length),"merged peaks\n")
})
