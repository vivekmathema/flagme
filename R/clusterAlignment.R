

setMethod("compress","clusterAlignment",
function(object,verbose=TRUE,...) {
  for(i in 1:length(object@alignments)) {
	pA<-object@alignments[[i]]
    object@alignments[[i]]<-compress(pA,verbose=FALSE)
  }
  new("clusterAlignment",object)
})

setMethod("decompress","clusterAlignment",
function(object,verbose=TRUE,...) {
  for(i in 1:length(object@alignments)) {
	pA<-object@alignments[[i]]
    object@alignments[[i]]<-decompress(pA,verbose=FALSE)
  }
  new("clusterAlignment",object)
})

clusterAlignment <- function(pD, runs=1:length(pD@rawdata),
                             timedf=NULL, usePeaks=TRUE, verbose=TRUE,...) {
    n <- length(runs)
    if(usePeaks)
        nr <- length(pD@peaksdata)
    else
        nr <- length(pD@rawdata)
    alignments <- vector("list", n*(n-1)/2)
    aligned <- matrix(-1, nr, nr)
    colnames(aligned) <- names(pD@rawdata)
    rownames(aligned) <- names(pD@rawdata)
    dist <- matrix(0, n, n)
    colnames(dist) <- names(pD@rawdata)[runs]
    rownames(dist) <- names(pD@rawdata)[runs]
    count <- 0
    ## browser()
    for(i in 1:(n-1))
    {
        run.i <- runs[i]
        for(j in (i+1):n)
        {
            run.j <- runs[j]
            count <- count+1
            if(verbose)
            {
                cat("[clusterAlignment] Aligning",
                    names(pD@rawdata)[run.i], "to",
                    names(pD@rawdata)[run.j], "\n")
            }
            if(usePeaks)
            {
                alignments[[count]] <-
                    peaksAlignment(pD@peaksdata[[run.i]],
                                   pD@peaksdata[[run.j]],
                                   pD@peaksrt[[run.i]],
                                   pD@peaksrt[[run.j]],
                                   usePeaks=usePeaks,
                                   timedf=timedf[[count]],
                                   verbose=verbose, ...)
             }
            else
            {
                alignments[[count]] <-
                    peaksAlignment(pD@rawdata[[run.i]],
                                   pD@rawdata[[run.j]],
                                   pD@rawrt[[run.i]],
                                   pD@rawrt[[run.j]],
                                   usePeaks=usePeaks, timedf=NULL,
                                   verbose=verbose, ...)
            }
            aligned[runs[i],runs[j]] <- aligned[runs[j],runs[i]] <- count
            dist[j,i] <- dist[i,j] <- alignments[[count]]@dist
	}
    }
    merge <- hclust(as.dist(dist), method="average")$merge
    merge.copy <- merge
    for(i in 1:length(runs))
        {merge[which(merge.copy == (-i))] <- (-runs[i])}
    new("clusterAlignment", runs=runs, aligned=aligned,
        gap=alignments[[1]]@gap, D=alignments[[1]]@D, dist=dist,
        alignments=alignments, merge=merge)
}

setMethod("show","clusterAlignment",
function(object) {
   cat("An object of class \"",class(object),"\"\n",sep="")
   cat("Pairwise distance matrix\n")
   print(object@dist)
   #cat(length(object@mz), "m/z bins - range: (",range(object@mz),")\n",sep=" ")
   #cat("scans:",sapply(object@rawdata,ncol),"\n",sep=" ")
   #cat("peaks:",sapply(object@peaksdata,ncol),"\n",sep=" ")
})

setMethod("plot","clusterAlignment",function(x,y,...) .plotcA(x,...))

.plotcA<-function(object,alignment=1,...) {
    rn<-rownames(object@aligned)
    for(i in alignment) {
	  ind<-which(object@aligned==i,arr.ind=TRUE)[2,]
	  plot(object@alignments[[i]],main=paste("D=",object@D," gap=",object@gap,sep=""),xlab=paste("Peaks ",rn[ind[1]],sep=" - "),ylab=paste("Peaks ",rn[ind[2]],sep=" - "),...)
	}
}

