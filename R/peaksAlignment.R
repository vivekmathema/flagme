setMethod("compress","peaksAlignment",
  function(object,verbose=TRUE,...) {
    if(object@compressed) {
      if (verbose)
  	    cat("[compress.peaksAlignment] Already compressed.\n")
      return(object)
    }
    #if (verbose)
    #  before<-sum(ll(object)$KB)
    object@r<-as.matrix.csc(1-object@r)
    #if (verbose) {
    #  after<-sum(ll(object)$KB)
  	#  cat("[compress.peaksAlignment] Object is",(before-after),"KB smaller.\n")
    #}
    object@compressed<-TRUE
    new("peaksAlignment",object)
  })

setMethod("decompress","peaksAlignment",
  function(object,verbose=TRUE,...) {
    if(!(object@compressed)) {
      if(verbose)
 	    cat("[decompress.peaksAlignment] Already decompressed.\n")
      return(object)
    }
    #if (verbose)
    #  before<-sum(ll(object)$KB)
    object@r<-1-(as.matrix(object@r))
    #if (verbose) {
    #  after<-sum(ll(object)$KB)
  	#  cat("[decompress.peaksAlignment] Object is",(after-before),"KB larger.\n")
    #}
    object@compressed<-FALSE
    new("peaksAlignment",object)
  })


#peaksAlignment<-function(d1,d2,t1,t2,D=1000,timedf=NULL,gap=.5,df=30,w=25,mzind=24,metric=1,transform="none",verbose=TRUE,usePeaks=TRUE,compress=TRUE) {
peaksAlignment<-function(d1,d2,t1,t2,gap=.5,D=1000,timedf=NULL,df=30,verbose=TRUE,usePeaks=TRUE,compress=TRUE) {

  #r<-switch(metric,
           #normDotProduct(d1,d2,t1,t2,D=D,df=df+abs(ncol(d1)-ncol(d2)),timedf=timedf,verbose=verbose),  # metric=1
           #cos.ndp.rank(d1,d2,t1,t2,D=D,timedf=timedf,verbose=verbose),                          # metric=2
           #euclidean(d1,d2,t1,t2,D=D,timedf=timedf,verbose=verbose),                             # metric=3
           #soai(d1,d2,t1,t2,D=D,timedf=timedf,verbose=verbose),                                  # metric=4
           #manhattan(d1,d2,t1,t2,D=D,timedf=timedf,verbose=verbose),                             # metric=5
           #window.metric.c(d1,d2,mzind=mzind,df=df,w=w))                                         # metric=6
		   #)
		   
  r<-normDotProduct(d1,d2,t1,t2,D=D,df=df+abs(ncol(d1)-ncol(d2)),timedf=timedf,verbose=verbose)

  if (verbose)
    cat("[peaksAlignment] Comparing",ncol(d1),"peaks to",ncol(d2),"peaks -- gap=", gap, "D=", D, "\n")
  v<-dp(r,gap=gap,verbose=verbose)
  sim<-0
  for(i in 1:nrow(v$match)) {
    sim<-sim+r[v$match[i,1],v$match[i,2]]
  }
  sim<-sim/nrow(v$match)
  if (verbose)
    cat("[peaksAlignment] ",nrow(v$match),"matched.  Similarity=",sim,"\n")
  object<-new("peaksAlignment",v=v,r=r,dist=sim,compressed=FALSE,gap=gap,D=D)
  if(compress) {
    compress(object)
  } else {
    object
  }
}


setMethod("show","peaksAlignment",
function(object) {
   cat("An object of class \"",class(object),"\"\n",sep="")
   dm<-dim(object@r)
   cat(dm[1],"peaks against",dm[2],"peaks:",nrow(object@v$match),"matched\n")
})


setMethod("plot","peaksAlignment",function(x,y,...) {require(gplots); .plotpA(x,...)})

.plotpA<-function(object,xlab="Peaks - run 1",ylab="Peaks - run 2",plotMatches=TRUE,matchPch=19,matchLwd=3,matchCex=.5,matchCol="black",col=colorpanel(50,"black","blue","white"),breaks=seq(0,1,length=51),...) {
		if(object@compressed)
		  object<-decompress(object,...)
	r<-object@r
	image(1:nrow(r),1:ncol(r),r,col=col,xlab=xlab,ylab=ylab,...)
	if (plotMatches)
	  points(object@v$match,pch=matchPch,col=matchCol,cex=matchCex,type="b",lwd=matchLwd)
}

