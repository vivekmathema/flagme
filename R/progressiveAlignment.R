
setMethod("decompress", "progressiveAlignment",
          function(object, verbose=TRUE, ...){
              if(object@merges[[1]]$compressed == FALSE) {
                  if(verbose)
                      cat("[decompress.progressiveAlignment] Already decompressed.\n")
                  return(object)
              }
              for(i in 1:length(object@merges)) {
                  if(object@merges[[i]]$compressed) {
                      object@merges[[i]]$r <- 1-as.matrix(object@merges[[i]]$r)
                      object@merges[[i]]$compressed <- FALSE
                  }
              }
              new("progressiveAlignment", object)
          })

setMethod("compress","progressiveAlignment",
          function(object,verbose=TRUE,...) {
              if(object@merges[[1]]$compressed) {
                  if(verbose)
                      cat("[compress.progressiveAlignment] Already compressed.\n")
                  return(object)
              }
              for(i in 1:length(object@merges)) {
                  if(!(object@merges[[i]]$compressed)) {
                      object@merges[[i]]$r <- as.matrix.csc(1-object@merges[[i]]$r)
                      object@merges[[i]]$compressed <- TRUE
                  }
              }
              new("progressiveAlignment",object)
          })


progressiveAlignment <- function(pD, cA, D=50, gap=0.5, verbose=TRUE,
                                 usePeaks=TRUE, df=30, compress=TRUE,
                                 type=2){
    m <- cA@merge
    merges <- vector("list", nrow(m))
    if(usePeaks)
        pd <- pD@peaksdata
    else
        pd <- pD@rawdata	
    pD <- NULL
    cA <- decompress(cA, verbose=verbose)
    
    for(i in 1:nrow(m)){
        if(verbose)
            cat("[progressiveAlignment] Doing merge", m[i,], "\n")
	if(m[i,1] < 0)
        {
            left.runs <- (-m[i,1])
            left.ind <- matrix(1:ncol(pd[[left.runs]]), ncol=1)
	}
        else
        {
            left.runs <- merges[[m[i,1]]]$runs
            left.ind <- merges[[m[i,1]]]$ind
	}
        
	if(m[i,2] < 0)
        {
            right.runs <- abs(m[i,2])
            right.ind <- matrix(1:ncol(pd[[right.runs]]), ncol=1)
	}
        else
        {
            right.ind <- merges[[m[i,2]]]$ind
            right.runs <- merges[[m[i,2]]]$runs
	}
	
	if(verbose){
            cat("[progressiveAlignment] left.runs:", left.runs, ", ")
            cat("right.runs:", right.runs, "\n")
	}
	
	if(length(right.runs) == 1 & length(left.runs) == 1)
        {
            al <- cA@alignments[[cA@aligned[left.runs,right.runs]]]
            v <- al@v
            sim <- al@r
            mi <- .merge.indices(nrow(left.ind), nrow(right.ind), v$match)
            new.ind <- mi
	}
        else
        {
            sim <- matrix(0, nrow=nrow(left.ind), ncol=nrow(right.ind))
            count <- matrix(0, nrow=nrow(left.ind), ncol=nrow(right.ind))
            if (verbose)
                cat("[progressiveAlignment] (dot=50) going to",
                    nrow(sim), ":")
            for(r in 1:nrow(sim)){
		if(verbose)
                    if(r %%50==0) cat(".")
		for(c in max(1, r-df) : min(r+df, ncol(sim))){
                    
                    for(lr in 1:length(left.runs)){
			for(rr in 1:length(right.runs)){
                            ai <- cA@aligned[left.runs[lr], right.runs[rr]]
                            #this.sim<-cA$alignments[[ai]]$r
                            
                            lind <- left.ind[r,lr]
                            rind <- right.ind[c,rr]
                            if(!is.na(lind) & !is.na(rind)){
                                if(left.runs[lr] < right.runs[rr])
                                    sim[r,c] <- sim[r,c] + cA@alignments[[ai]]@r[lind,rind]
				else
                                    sim[r,c] <- sim[r,c] + cA@alignments[[ai]]@r[rind,lind]
				count[r,c] = count[r,c]+1
                            }			  
                            
			}
                    }
                    
		}
            }
            if(verbose)
                cat("\n")
            
            if(type == 2) # RR
            {
                v <- dynRT(S=sim)
                v$match <- v$match[!is.na(v$match[,2]),] # remove non-matched peaks
            }
            if(type == 1)
            {
                sim[count > 0] <- sim[count>0]/count[count>0]
                sim[sim == 0] <- 1
                v <- dp(sim, gap=gap, verbose=verbose)#
            }    
            mi <- .merge.indices(nrow(left.ind), nrow(right.ind), v$match)
            new.ind <- cbind(left.ind[mi[,1],],right.ind[mi[,2],])
	}
	rownames(new.ind) <- 1:nrow(new.ind)
	merges[[i]] <- list(ind=new.ind, mi=mi,
                            runs=c(left.runs,right.runs),
                            left=left.runs, right=right.runs, r=sim, compressed=FALSE)
    }
    cA <- NULL
    if(verbose)
        print(gc())
    pA <- new("progressiveAlignment", merges=merges)
    if(compress)
        return(compress(pA, verbose=verbose))
    else
	return(pA)
}


.merge.indices <- function(nl, nr, m) {
    lind <- cbind(1:nl, 0)
    lind[lind[,1] %in% m[,1],2] <- 1
    rind <- cbind(1:nr,0)
    rind[rind[,1] %in% m[,2],2] <- 1
    mg <- matrix(NA, nrow=nrow(lind)+nrow(rind)-nrow(m), ncol=2)
    #print(dim(mg))
    li <- 1; ri <- 1; i <- 1
    while(i <= nrow(mg)){
        if(li > nrow(lind) & ri <= nrow(rind))
        {
            mg[i,] <- c(NA,rind[ri,1])
            ri <- ri+1
            i <- i+1
            next
	}
        if(li <= nrow(lind) & ri > nrow(rind))
        {
            mg[i,] <- c(lind[li,1], NA)
            li <- li+1	
            i <- i+1    
            next
	}
        if(lind[li,2] == 1)
        {
            if(rind[ri,2] == 1)
            {  # match
                mg[i,] <- c(lind[li,1], rind[ri,1])
		ri <- ri+1
		li <- li+1
		#cat("match",i,li,ri,"-",mg[i,],"\n")
            }
            else
            {             # right unmatched
                mg[i,] <- c(NA, rind[ri,1])
		ri <- ri+1
		#cat("right unmatched",i,li,ri,"-",mg[i,],"\n")
            }
	}
        else
        {
            if(rind[ri,2] == 1)
            {  # left unmatched
                mg[i,] <- c(lind[li,1], NA)
		li <- li+1	    
		#cat("left unmatched",i,li,ri,"-",mg[i,],"\n")
            }
            else
            {             # both unmatched
                mg[i,] <- c(lind[li,1], NA)
		#cat("both unmatched - A",i,li,ri,"-",mg[i,],"\n")
		i <- i+1
                mg[i,] <- c(NA, rind[ri,1])
		li <- li+1
		ri <- ri+1
		#cat("both unmatched - B",i,li,ri,"-",mg[i,],"\n")
            }
	}
	i <- i+1
    }
    mg
}

setMethod("show","progressiveAlignment",
function(object){
   cat("An object of class \"", class(object), "\"\n", sep="")
   cat(length(object@merges), "merges\n")
})


