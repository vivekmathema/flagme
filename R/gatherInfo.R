#' Gathers abundance informations from an alignment
#' 
#' Given an alignment table (indices of matched peaks across several samples)
#' such as that within a \code{progressiveAlignment} or
#' \code{multipleAlignment} object, this routines goes through the raw data and
#' collects the abundance of each fragment peak, as well as the retention times
#' across the samples.
#' 
#' This procedure loops through the the table of matched peaks and gathers the
#' 
#' @param pD a \code{peaksDataset} object, to get the abundance data from
#' @param obj either a \code{multipleAlignment} or \code{progressiveAlignment}
#' object
#' @param newind list giving the
#' @param method method used to gather abundance information, only \code{apex}
#' implemented currently.
#' @param findmzind logical, whether to take a subset of all m/z indices
#' @param useTIC logical, whether to use total ion current for abundance
#' summaries
#' @param top only use the top \code{top} peaks
#' @param intensity.cut percentage of the maximum intensity
#' @return Returns a list (of lists) for each row in the alignment table.  Each
#' list has 3 elements: \item{mz}{a numerical vector of the m/z fragments used}
#' \item{rt}{a numerical vector for the exact retention time of each peak
#' across all samples} \item{data}{matrix of fragment intensities.  If
#' \code{useTIC = TRUE}, this matrix will have a single row}
#' @author Mark Robinson
#' @seealso \code{\link{imputePeaks}}
#' @references Mark D Robinson (2008).  Methods for the analysis of gas
#' chromatography - mass spectrometry data \emph{PhD dissertation} University
#' of Melbourne.
#' @keywords manip
#' @examples
#' 
#'   require(gcspikelite)
#' 
#'   ## paths and files
#'   gcmsPath <- paste(find.package("gcspikelite"), "data", sep = "/")
#'   cdfFiles <- dir(gcmsPath, "CDF", full = TRUE)
#'   eluFiles <- dir(gcmsPath, "ELU", full = TRUE)
#' 
#'   ## read data, peak detection results
#'   pd <- peaksDataset(cdfFiles[1:2], mz = seq(50, 550), rtrange = c(7.5, 8.5))
#'   pd <- addAMDISPeaks(pd, eluFiles[1:2])
#' 
#'   ## multiple alignment
#'   ma <- multipleAlignment(pd, c(1,1), wn.gap = 0.5, wn.D = 0.05, bw.gap = 0.6, 
#'                           bw.D = 0.2, usePeaks = TRUE, filterMin = 1, df = 50,
#'                           verbose = TRUE, metric = 1, type = 1)
#' 
#'   ## gather apex intensities
#'   d <- gatherInfo(pd, ma)
#' 
#'   ## table of retention times
#'   nm <- list(paste("MP", 1:length(d), sep = ""), c("S1", "S2"))
#'   rts <- matrix(unlist(sapply(d, .subset, "rt")), byrow = TRUE, nc = 2, 
#'                 dimnames = nm)
#' 
#' @export gatherInfo
gatherInfo <- function(pD, obj, newind=NULL, method=c("apex"),
                       findmzind=TRUE, useTIC=FALSE, top=NULL,
                       intensity.cut=.05){
    ## ## correlationAlignment() data implementation
    ## ## start 
    ## if(class(obj) == 'correlationAlignment')
    ## {
    ##     ## slot rt
    ##     l.rt <- lapply(1:nrow(obj@Alignment), function(z){
    ##         sapply(seq(along=obj@Alignment), function(x, y){
    ##             n <- names(obj@Alignment)[x] #!
    ##             r <- pD@peaksrt[[n]][obj@Alignment[y,x]]
    ##         }, y=1:nrow(obj@Alignment))[z,]
    ##     })
    ##     names(l.rt) <- sapply(l.rt, function(x){names(x) <- 'rt'})    
    ##     ## slot data
    ##     l.d <- lapply(1:nrow(obj@Alignment), function(y){
    ##         rx <- lapply(1:ncol(obj@Alignment), function(x){
    ##             nam <- names(obj@Alignment)[x]
    ##             idx <- obj@Alignment[y,x]
    ##             if(is.na(idx) == FALSE){
    ##                 ## dat <- pD@peaksdata[nam][1][[1]][,idx]
    ##                 ## must match
    ##                 ## the column name not the clolumn position
    ##                 dat <- pD@peaksdata[nam][1][[1]][,which(colnames(pD@peaksdata[nam][1][[1]]) == idx)]
    ##             }else{
    ##                 dat <- rep(NA, length(pD@peaksdata[nam][1][[1]][,idx]))
    ##             }
    ##             return(dat)
    ##         })
    ##         dd <- do.call(cbind, rx)
    ##         return(dd)
    ##     })
    ##     names(l.d) <- sapply(l.d, function(x){names(x) <- 'data'})
    ##     ## trace back raw id spectra
    ##     l.Raw <- lapply(1:nrow(obj@Alignment), function(x){
    ##         cbind.data.frame(Sample=colnames(obj@Alignment),
    ##                          rawID=as.numeric(obj@Alignment[x,]))
    ##     })
    ##     names(l.Raw) <- sapply(l.Raw, function(x){names(x) <- 'traceRaw'})
        
    ##     res <- lapply(1:length(l.d), function(x){
    ##         c(l.Raw[x], l.rt[x], l.d[x]) 
    ##     })
    ##     names(res) <- paste("feat", seq(along=res), sep='_')
    ##     return(res)
    ## }
    ##     ## end

    ## pind - indices of peaks
    ## rind - indices of runs
    ## newind ... comes from an imputation step, if at all.
    method <- match.arg(method)
    if(is(obj, "multipleAlignment"))
    {
        pind <- obj@betweenAlignment@ind
        rind <- obj@betweenAlignment@runs
        groups <- obj@betweenAlignment@groups
    }
    else
        if(is(obj, "progressiveAlignment"))
        {
            n <- length(obj@merges)
            pind <- obj@merges[[n]]$ind   # for example, if obj=ma$pas[1]
            rind <- obj@merges[[n]]$runs
            groups <- rep(names(obj),length(rind))
        }

    ## collect rt
    out <- matrix(NA, nrow(pind), ncol(pind))
    for(j in 1:ncol(pind)){
        data <- pD@peaksrt[[rind[j]]]
        keep <- which(!is.na(pind[,j]))
        out[keep,j] <- data[pind[keep,j]]
    }
    colnames(out) <- paste(groups,rind,sep=".")
    rownames(out) <- 1:nrow(out)

    ## collect mz and intensity
    peaks <- vector("list", nrow(pind))
    for(i in 1:length(peaks)){
        if(findmzind)
        {
            mzind <- rep(0, length(pD@mz))
            for(j in 1:ncol(pind)){
                data <- pD@peaksdata[[rind[j]]]
                if(!is.na(pind[i,j]))
                    mzind <- mzind + data[,pind[i,j]]
            }
            ## mzind <- which(mzind > 0)
            mzind <- which(mzind > intensity.cut * max(mzind)) # intensity filter
        }
        else
        {
            mzind <- 1:length(pD@mz)
        }
        if(useTIC)
            this.d <- matrix(, nrow=1, ncol=ncol(pind))
        else
            this.d <- matrix(, nrow=length(mzind), ncol=ncol(pind))
        
        for(j in 1:ncol(pind)){
            r <- rind[j]
            overlap <- NULL
            if(!is.na(pind[i,j]))
                if(method == "apex")
                {
                    if(useTIC)
                    {
                        this.d[1,j] <- sum(pD@rawdata[[r]][,pD@peaksind[[r]][pind[i,j]]])
                    } ##
                    else
                    {
                        ## this.d[,j] <- pD@rawdata[[r]][mzind,pD@peaksind[[r]][pind[i,j]]] # ?
                        this.d[,j] <- pD@peaksdata[[r]][mzind,pind[i,j]]
                    }
                }
                else
            {
                ## method = "area"
                ## nn <- length(pD@peaksind.start[[r]])
                ## w <- which(round(pD@peaksind.start[[r]][-1]) < round(pD@peaksind.end[[r]][-nn]))
                ## overlap <- cbind(w,w+1)
                ## v <- which(pind[i,j] == overlap, arr.ind=TRUE)
                ## if(nrow(v) > 0) {
                ##     ## do mixture thing. 
                ##     pks <- overlap[v[1,1],]
                ##     cat("overlap -- run",r,"-- peaks",pks,"--",v[1,2],"\n")
                ##     if(r == 12 & pks[1] == 4)
                ##         return(list(mzind=mzind))
                ##     this.d[,j] <- partitionWithNormalMixture2(pD,r,pks[1],pks[2],mzind,useTIC=useTIC)[,v[1,2]]
                ##     if(j == ncol(pind)) cat("\n")
                ## }else{
                ##     ## just simple sum of intensity
                ##     st <- pD@peaksind.start[[r]][pind[i,j]]
                ##     en <- pD@peaksind.end[[r]][pind[i,j]]
                ##     if (useTIC)
                ##         this.d[1,j] <- sum(pD@rawdata[[r]][,st:en])  ##
                ##     else
                ##         this.d[,j] <- rowSums(pD@rawdata[[r]][mzind,st:en])
                ## }
            }
            if(!is.null(newind))
                if(!is.na(newind$apex[i,j])){
                    if(method == "apex"){
                        if(useTIC)
                            this.d[1,j] <- sum(pD@rawdata[[r]][,newind$apex[i,j]])  ##
                        else
                            this.d[,j] <- pD@rawdata[[r]][mzind,newind$apex[i,j]]
                        out[i,j] <- pD@rawrt[[r]][newind$apex[i,j]]
                    }else{ ## method="area"
                        ## st <- newind$start[i,j]
                        ## en <- newind$end[i,j]
                        ## cat(r,i,j,st,en,"\n")
                        ## if (useTIC)
                        ##     this.d[1,j] <- sum(pD@rawdata[[r]][,st:en])  ##
                        ## else
                        ##     this.d[,j] <- rowSums(pD@rawdata[[r]][mzind,st:en])
                    }
                }
        }
        colnames(this.d) <- paste(groups, rind, sep=".")
        peaks[[i]] <- list(rt=out[i,], mz=pD@mz[mzind], data=this.d)
    }
    if(!is.null(top)){
        for(i in seq(along=peaks)){
            o <- sort(order(-rowSums(peaks[[i]]@data, na.rm=TRUE))[1:top])
            peaks[[i]]@mz <- peaks[[i]]@mz[o]
            peaks[[i]]@data <- peaks[[i]]@data[o,]
        }
    }
    return(peaks)
}

