
## peaksDataset<-function(fns=dir(,"[Cc][Dd][Ff]"),verbose=TRUE,mz=seq(50,550),rtDivide=60,rtrange=NULL) {
##   rawdata<-vector("list",length(fns))
##   rawrt<-vector("list",length(fns))
##   for(i in 1:length(fns)) {
##     if (verbose)
##       cat(" Reading ", fns[i],"\n")
##     a<-xcmsRaw(fns[i])
## 	if (is.null(mz))
## 	  mz<-seq(a@mzrange[1],a@mzrange[2])
##     this.mz<-seq(a@mzrange[1],a@mzrange[2])
## 	rawrt[[i]]<-a@scantime/rtDivide
## 	nc<-length(rawrt[[i]])
##     if(length(rtrange)==2) {
##       w<-which(rawrt[[i]] >= rtrange[1] & rawrt[[i]] <= rtrange[2])
## 	  rawrt[[i]]<-rawrt[[i]][w]
## 	} else {
## 	  w<-seq(1,nc)
## 	}
## 	rawdata[[i]]<-a@env$profile[,w]
## 	rm(a)
## 	if (length(this.mz)!=length(mz)) {
## 	  d<-matrix(0,nrow=length(mz),ncol=length(w))
## 	  d[match(this.mz,mz),]<-rawdata[[i]]
## 	  rawdata[[i]]<-d
## 	}
##   }
##   gc()
##   nm<-lapply(fns,FUN=function(u) { sp<-strsplit(u,split="/")[[1]]; sp[length(sp)]})
##   nm<-sub(".CDF","",nm)
##   names(rawdata)<-names(rawrt)<-nm
##   #d<-list(mz=mz,files=fns,rawdata=rawdata,rawrt=rawrt)
##   #new("peaksDataset",d)
##   new("peaksDataset",rawdata=rawdata,rawrt=rawrt,mz=mz,files=fns)
## }

peaksDataset <- function(fns=dir(,"[Cc][Dd][Ff]"), verbose=TRUE,
                         mz=seq(50,550), rtDivide=60, rtrange=NULL)
{
    rawdata <- vector("list", length(fns))
    rawrt <- vector("list", length(fns))
    for(i in 1:length(fns))
    {
        if(verbose)
            cat(" Reading ", fns[i],"\n")
        a <- xcmsRaw(fns[i])
        if(is.null(mz) == TRUE)
        {
            mz <- seq(a@mzrange[1], a@mzrange[2])
            this.mz <- seq(a@mzrange[1], a@mzrange[2])
        }
        else 
        {
            this.mz <- mz
        }
        rawrt[[i]] <- a@scantime/rtDivide
        nc <- length(rawrt[[i]])
        if(length(rtrange) == 2)
        {
            w <- which(rawrt[[i]] >= rtrange[1] & rawrt[[i]] <= rtrange[2])
            rawrt[[i]] <- rawrt[[i]][w]
        }
        else
        {
            w <- seq(1, nc)
        }
        rawdata[[i]] <- a@env$profile[,w]
        rm(a) ## alleggerisco output
        if(length(this.mz)!= length(mz))
        {
            d <- matrix(0, nrow=length(mz), ncol=length(w))
            ## bug feb2015
            # il problema si verifica quando utilizzo un range di massa
            # inferiore a quello dell'acquisizione effettuata dallo
            # spettrometro. 
            #
            # provo ad escludere NA 
            m <- match(this.mz, mz)
            mcl <- m[which(m != is.na(match(this.mz, mz)))]
            d[mcl,] <- rawdata[[i]] ## bug
            ###
            ## d[match(this.mz, mz),] <- rawdata[[i]] ## bug
            rawdata[[i]] <- d
        }
    }
    gc()
    nm <- lapply(fns, FUN=function(u){
        sp <- strsplit(u, split="/")[[1]]; sp[length(sp)]
    }
                 )
    nm <- sub(".CDF", "", nm)
    names(rawdata) <- names(rawrt) <- nm
    # d <- list(mz=mz, files=fns, rawdata=rawdata, rawrt=rawrt)
    # new("peaksDataset",d)
    new("peaksDataset", rawdata=rawdata, rawrt=rawrt, mz=mz, files=fns)
}


.plotpD<-function(object,runs=1:length(object@rawdata),mzind=1:nrow(object@rawdata[[1]]),mind=NULL,plotSampleLabels=TRUE,calcGlobalMax=FALSE,peakCex=0.8,plotPeaks=TRUE,plotPeakBoundaries=FALSE,plotPeakLabels=FALSE,plotMergedPeakLabels=TRUE,mlwd=3,usePeaks=TRUE,plotAcrossRuns=FALSE,overlap=F,rtrange=NULL,cols=NULL,thin=1,max.near=median(object@rawrt[[1]]),how.near=50,scale.up=1,...) {
  # only specific one of ma or runs
  if(is.null(rtrange))
    rtrange<-c(min(object@rawrt[[1]]),max(object@rawrt[[1]]))
  if(is.factor(cols))
    cols<-c("black","green","blue","red","grey","orange")[cols]
  if(is.null(cols))
    cols<-rep(c("black","green","blue","red"),each=length(object@rawdata)/4)
  if(length(cols)!=length(runs))
    cols<-rep("black",length(runs))
  if (usePeaks) {
    datalist<-object@peaksdata
	rtlist<-object@peaksrt
  } else {
    datalist<-object@rawdata
	rtlist<-object@rawrt
  }
  if(calcGlobalMax) {
    mx<-0
    for(i in 1:length(runs)) {
      ind<-runs[i]
	  mxind<-which(abs(object@rawrt[[ind]]-max.near)<how.near)
  	  b<-colSums(matrix(object@rawdata[[ind]][mzind,],nrow=length(mzind)))
	  if (mx < max(b[mxind]))
	    mx<-max(b[mxind])
	}
  }
  if( length(object@peaksdata)!=length(object@rawdata) )
    plotPeaks<-FALSE
  for(i in 1:length(runs)) {
    ind<-runs[i]
	b<-colSums(matrix(object@rawdata[[ind]][mzind,],nrow=length(mzind)))
	if (!calcGlobalMax) {
	  mxind<-which(abs(object@rawrt[[ind]]-max.near)<how.near)
	  mx<-max(b[mxind])
	}
    b<-b/mx*scale.up
	if (!overlap) {
      if(i==1) {
	    plot(object@rawrt[[ind]],b,type="l",xlab="Retention Time",ylab="",xlim=rtrange,ylim=c(-.5,length(runs)),col=cols[i],xaxs="i",axes=FALSE,...)
            axis(1)
	  } else {
	    lines(object@rawrt[[ind]],b+i-1,col=cols[i],...)
	  }
	  if (plotSampleLabels)
	    text(rtrange[1]+.1*diff(rtrange),i-.5,names(object@rawdata)[ind],col="black")
    } else {
      if(i==1) {
	    plot(object@rawrt[[ind]],b,type="l",xlab="Retention Time",ylab="",xlim=rtrange,ylim=c(0,1),col=cols[i],xaxs="i",axes=FALSE)
            axis(1)
	  } else {
	    lines(object@rawrt[[ind]],b,col=cols[i])
	  }
	}
    if(plotPeaks & !overlap) {
   	  rt<-rtlist[[ind]]; d<-min(diff(rt))*.4
	  for(j in 1:length(rt)) {
	    if(rt[j] > rtrange[1] & rt[j] < rtrange[2]) {
		  if (plotPeakBoundaries) {
			st<-object@rawrt[[ind]][ object@peaksind.start[[ind]][j] ]
			en<-object@rawrt[[ind]][ object@peaksind.end[[ind]][j] ]
                           #cat(ind,st,en,rt[j],"\n")
			if(j%%2==0)
	          lines(c(st,st,rt[j],rt[j],rt[j],en,en),c(i-1+.025,i-1,i-1,i-1+.3,i-1,i-1,i-1+.025),col="blue")
			else
	          lines(c(st,st,rt[j],rt[j],rt[j],en,en),c(i-1+.025,i-1,i-1,i-1+.3,i-1,i-1,i-1+.025)-.025,col="red")
	        #lines(c(rt[j],rt[j]),c(i-1,i-1+.3))
		  } else {
	        lines(c(rt[j]-d,rt[j]+d),c(i-1,i-1))
	        lines(c(rt[j],rt[j]),c(i-1,i-1+.3))
		  }
	    }
	  }
	  if(plotPeakLabels)
	    text(rt,i-1+.15,labels=1:length(rt),adj=0,cex=peakCex)
    }
  }
  if(!is.null(mind) & !overlap) {
	  rws<-as.integer(seq(1,nrow(mind),by=thin))
	  for(i in rws) {
	    if( !is.na(mind[i,1]) & plotMergedPeakLabels )
		  text(rtlist[[runs[1]]][mind[i,1]],0,srt=90,i,col="blue",pos=1,cex=peakCex*.8,adj=1)
	    for(j in 1:(ncol(mind)-1)) {
	      ind1<-mind[i,j]
	      ind2<-mind[i,j+1]
		  if(!is.na(ind1) & !is.na(ind2)) {
		    rt1<-rtlist[[runs[j]]][ind1]
			rt2<-rtlist[[runs[j+1]]][ind2]
			#cat(ind1,ind2,runs[j],runs[j+1],rt1,rtrange[1],rt2,rtrange[2],"\n")
			if( (rt1 > rtrange[1] | rt2 > rtrange[1]) & (rt1 < rtrange[2] | rt2 < rtrange[2])) {
			  if(j > 1)
		        lines(c(rt1,rt2),c(j+.4-1,j+.95-1),lwd=mlwd,lty=1,col="darkgrey")
			  else
		        lines(c(rt1,rt2),c(j+.4-1,j+.95-1),lwd=mlwd,lty=1,col="darkgrey")
			}
		  }
		}
		if (plotAcrossRuns) {
 		  if (ncol(mind) <= 2) next
	      for(j in 1:(ncol(mind)-2)) {
	        ind1<-mind[i,j]
	        ind2<-mind[i,j+1]
	        ind3<-mind[i,j+2]
		    if(!is.na(ind1) & is.na(ind2) & !is.na(ind3)) {
		      rt1<-rtlist[[runs[j]]][ind1]
			  rt3<-rtlist[[runs[j+2]]][ind3]
			  #cat(rt1,rtrange[1],rt3,rtrange[2],"\n")
			  if( (rt1 > rtrange[1] | rt3 > rtrange[1]) & (rt1 < rtrange[2] | rt3 < rtrange[2]))
		        lines(c(rt1,rt3),c(j+.4-1,j+.95),lwd=mlwd,lty=1,col="green")
		    }
		  }
		  if (ncol(mind) <= 3) next
	      for(j in 1:(ncol(mind)-3)) {
	        ind1<-mind[i,j]
	        ind2<-mind[i,j+1]
	        ind3<-mind[i,j+2]
	        ind4<-mind[i,j+3]
		    if(!is.na(ind1) & is.na(ind2) & is.na(ind3) & !is.na(ind4)) {
		      rt1<-rtlist[[runs[j]]][ind1]
			  rt4<-rtlist[[runs[j+3]]][ind4]
			  #cat(rt1,rtrange[1],rt4,rtrange[2],"\n")
			  if( (rt1 > rtrange[1] | rt4 > rtrange[1]) & (rt1 < rtrange[2] | rt4 < rtrange[2]))
		        lines(c(rt1,rt4),c(j+.4-1,j+1.95),lwd=mlwd,lty=1,col="red")
		    }
		  }
	    }
	  }
  }
}

#setGeneric("plot", function(object, ...) standardGeneric("plot"))
#setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#setMethod("plot",signature("peaksDataset","missing"),function(x,y,...) plotpD(x,...))
setMethod("plot","peaksDataset",function(x,y,...) .plotpD(x,...))

setMethod("show","peaksDataset",
function(object) {
        cat("An object of class \"",class(object),"\"\n",sep="")
        cat(length(object@rawdata), "samples:",names(object@rawdata),"\n",sep=" ")
        cat(length(object@mz), "m/z bins - range: (",range(object@mz),")\n",sep=" ")
        cat("scans:",sapply(object@rawdata,ncol),"\n",sep=" ")
        cat("peaks:",sapply(object@peaksdata,ncol),"\n",sep=" ")
}) 


addAMDISPeaks<-function(object,fns=dir(,"[Eu][Ll][Uu]"),verbose=TRUE,...) {
  if(length(fns)!=length(object@rawdata))
    stop("Number of files must be the same as the number of runs (and must match).")
  mn<-1e6; mx<--1
  for(i in 1:length(object@rawrt)) {
    rng<-range(object@rawrt[[i]])
  	if(mn > rng[1]) mn<-rng[1]
	if(mx < rng[2]) mx<-rng[2]
  }
  if (verbose)
    cat("Reading retention time range:",c(mn,mx),"\n")
  for(i in 1:length(object@files)) {
    if (verbose)
	  cat("Reading",fns[i],"...")
    csin<-parseELU(fns[i],mz=object@mz,rtrange=c(mn,mx),...)
	thisind<-rep(NA,nrow(csin$tab))
	thisind.st<-rep(NA,nrow(csin$tab))
	thisind.en<-rep(NA,nrow(csin$tab))
	thisrt<-rep(NA,nrow(csin$tab))
	thisdata<-matrix(0,nrow=length(object@mz),ncol=nrow(csin$tab))
	for(j in 1:length(thisind)) {
	  thisind[j]<-which.min(abs(object@rawrt[[i]]-csin$tab[j,2]))
	  thisind.st[j]<-csin$tab[j,3]
	  thisind.en[j]<-csin$tab[j,4]
	  thisrt[j]<-object@rawrt[[i]][thisind[j]]
	  take<-which(csin$peaks[,j]>0)
	  #thisdata[take,j]<-sqrt(object@rawdata[[i]][take,thisind[j]])
	  thisdata[take,j]<-object@rawdata[[i]][take,thisind[j]]
	}
    if (verbose)
	  cat(" Done.\n")
	object@peaksind[[i]]<-thisind
	thisind.st[thisind.st==0]<-1
	object@peaksind.start[[i]]<-thisind.st
	object@peaksind.end[[i]]<-thisind.en
	object@peaksrt[[i]]<-thisrt
	object@peaksdata[[i]]<-thisdata
  }
  names(object@peaksdata)<-names(object@peaksrt)<-names(object@peaksind)<-names(object@rawdata)
  new("peaksDataset",object)
}

addChromaTOFPeaks<-function(object,fns=dir(,"[Tt][Xx][Tx]"),rtDivide=60,verbose=TRUE,...) {
  if(length(fns)!=length(object@rawdata))
    stop("Number of files must be the same as the number of runs (and must match).")
  mn<-1e6; mx<--1
  for(i in 1:length(object@rawrt)) {
    rng<-range(object@rawrt[[i]])
  	if(mn > rng[1]) mn<-rng[1]
	if(mx < rng[2]) mx<-rng[2]
  }
  if (verbose)
    cat("Reading retention time range:",c(mn,mx),"\n")
  for(i in 1:length(object@files)) {
    if (verbose)
	  cat("Reading",fns[i],"...")
    csin<-parseChromaTOF(fns[i],mz=object@mz,rtrange=c(mn,mx),rtDivide=rtDivide,...)
	thisrt<-csin$tab$rt
	thisind<-rep(NA,length(thisrt))
	thisdata<-matrix(0,nrow=length(object@mz),ncol=nrow(csin$tab))
	for(j in 1:length(thisind)) {
	  thisind[j]<-which.min(abs(object@rawrt[[i]]-csin$tab$rt[j]/rtDivide))
	  thisrt[j]<-object@rawrt[[i]][thisind[j]]
	  take<-which(csin$peaks[,j]>0)
	  #thisdata[take,j]<-sqrt(object@rawdata[[i]][take,thisind[j]])
	  thisdata[take,j]<-object@rawdata[[i]][take,thisind[j]]
	}
    if (verbose)
	  cat(" Done.\n")
	object@peaksind[[i]]<-thisind
	object@peaksrt[[i]]<-thisrt
	object@peaksdata[[i]]<-thisdata
  }
  names(object@peaksdata)<-names(object@peaksrt)<-names(object@peaksind)<-names(object@rawdata)
  new("peaksDataset",object)
}



setMethod("plotImage","peaksDataset",
  function(object,run=1,rtrange=c(11,13),main=NULL,mzrange=c(50,200),SCALE=log2,...) {
  mz<-object@mz
  for(i in run) {
    rt<-object@rawrt[[i]]
    d<-object@rawdata[[i]]	
    cols<-colorpanel(50,"green","blue","red")
    rind<-which(rt>rtrange[1]&rt<rtrange[2])
    mind<-which(mz>mzrange[1]&mz<mzrange[2])
    ff<-t(d[mind,rind])
    ff[ff==0]<-min(ff[ff>0])
	if( is.null(main) )
	  tit<-object@files[i]
	else
	  tit<-main
	#cat("min=",min(SCALE(ff),na.rm=TRUE)," max=",max(SCALE(ff),na.rm=T),"\n")
    image(rt[rind],mz[mind],SCALE(ff),col=cols,xlab="retention time",ylab="m/z",xlim=rtrange,main=tit,...)
    #abline(h=c(72.5,73.5,146.5,147.5),col="white",lty=3)
    #if (!is.null(pts)) {
    #  points(pts[,1],pts[,2],pch=19)
    #  matplot(t(pts[,8:9]),t(pts[,c(2,2)]),type="l",col="black",lty=1,add=TRUE)
    #}
  }
  return(invisible(ff))
})

