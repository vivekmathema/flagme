multipleAlignment <- function(pd, group, bw.gap = 0.8, wn.gap = 0.6,
                              bw.D = 0.20, wn.D = 0.05, filterMin = 1,
                              lite = FALSE, usePeaks = TRUE, df = 50,
                              verbose = TRUE, timeAdjust = FALSE,
                              doImpute = FALSE, metric = 2, type = 2, 
                              penality = 0.2){
  groups <- unique(group)
  cAs <- vector("list", length(groups))
  pAs <- vector("list", length(groups))
  timedf <- vector("list", length(groups))
  impute <- vector("list", length(groups))	

  ## do within alignment first
  for(i in 1:length(groups))
  {
    runs <- which(group == groups[i])
    if(timeAdjust == TRUE)
    {
      fullca <- clusterAlignment(pd, runs, gap = 0.05, D = 10, usePeaks = FALSE,
                                df = 100, verbose = verbose, metric = metric,
                                type = type)
      timedf[[i]] <- calcTimeDiffs(pd, fullca, verbose = verbose)
    }
      ## cAs[[i]] <- clusterAlignment(pd, runs, gap = wn.gap, D = wn.D,
      ##                              timedf = timedf[[i]], df = df,
      ##                              usePeaks = usePeaks,
      ##                              verbose = verbose, metric = 1, type = 1)  
      ## pAs[[i]] <- progressiveAlignment(pd, cAs[[i]], gap = wn.gap, D = wn.D,
      ##                                  df = df, usePeaks = usePeaks,
      ##                                  verbose = verbose, type = 1)
    cAs[[i]] <- clusterAlignment(pd, runs, gap = wn.gap, D = wn.D,
                                 timedf = timedf[[i]], df = df, 
                                 metric = metric, type = type, compress = TRUE,
                                 penality = penality) 
    pAs[[i]] <- progressiveAlignment(pd, cAs[[i]], gap = wn.gap, D = wn.D, 
                                     df = df, usePeaks = usePeaks, 
                                     compress = TRUE, type = type)
    if(doImpute)
    {
      if(timeAdjust)
      {
        impute[[i]] <- imputePeaks(pd, pAs[[i]], fullca, typ = 2)
      }
      else
      {
        impute[[i]] <- imputePeaks(pd, pAs[[i]], typ = 1)
      }
    }
  }

  names(cAs) <- groups
  names(pAs) <- groups

  if(length(groups) == 1)
  {
    pn <- length(pAs[[1]]@merges)
    ## fix this ... need a dummy object here
    wRA <- new("betweenAlignment", mergedPeaksDataset = new("peaksDataset"),
               ind = pAs[[1]]@merges[[pn]]$ind, runs = pAs[[1]]@merges[[pn]]$runs,
               imputeind = list(), filtind = list(), newind = list(),
               cA = cAs[[1]], pA = pAs[[1]],
               groups = as.character(rep(groups,nrow(cAs[[1]]@dist)))
               )
  } 
  else
  {
    wRA <- betweenAlignment(pd, cAs, pAs, impute, gap = bw.gap, D = bw.D,
                            filterMin = filterMin, usePeaks = usePeaks,
                            df = df, verbose = verbose, metric = metric,  
                            type = type, penality = penality)
  }
  
  ma <- new("multipleAlignment", clusterAlignmentList = cAs,
            progressiveAlignmentList = pAs, timeDiff = timedf,
            impute = impute, betweenAlignment = wRA)
  if(lite)
  {
    ma@clusterAlignmentList <- NULL
    ma@progressiveAlignmentList <- NULL
  }
  ma
}

setMethod("show","multipleAlignment",
          function(object){
    cat("An object of class \"", class(object), "\"\n", sep = "")
    cat(length(object@clusterAlignmentList), "groups:",
        sapply(object@progressiveAlignmentList, FUN = function(u)
          length(u@merges)+1), "samples, respectively.\n", sep = " ")
    cat(nrow(object@betweenAlignment@ind), "merged peaks","\n")
    ## cat(length(object@mz), "m/z bins - range: (",range(object@mz),")\n",sep=" ")
    ## cat("scans:",sapply(object@rawdata,ncol),"\n",sep=" ")
    ## cat("peaks:",sapply(object@peaksdata,ncol),"\n",sep=" ")
}
)

imputePeaks <- function(pD, obj, typ = 1, obj2 = NULL, filterMin = 1, 
                        verbose = TRUE){
  if(is(obj, "multipleAlignment")) 
  {
    ind <- obj@betweenAlignment@ind
    runs <- obj@betweenAlignment@runs
  } 
  else if(is(obj, "progressiveAlignment"))
  {
    n <- length(obj@merges)
    ind <- obj@merges[[n]]$ind
    runs <- obj@merges[[n]]$runs
  }
  newind <- matrix(,nrow = nrow(ind), ncol = ncol(ind))
  mask <- !is.na(ind) + 0
  rws <- which(rowSums(mask) >= filterMin & rowSums(mask) < ncol(mask))
  
  ind.scan <- ind
  for(i in 1:nrow(ind.scan)){  #  replace indices with scan numbers, then impute these 
    for(j in 1:ncol(ind.scan)){
      if(!is.na(ind.scan[i,j])) 
      {
        ind.scan[i,j] <- pD@peaksind[[runs[j]]][ind[i,j]]
      }
    } 
  }

  n <- ncol(ind.scan)
  approxs <- vector("list", n * (n - 1))
  count <- 0
  xy <- matrix(,n,n); rownames(xy) <- colnames(xy) <- 1:n
  
  if (typ == 1){  # impute from scan numbers
  
    for(i in 1:n){ 
      for(j in 1:n){ 
        if(i == j) next
        count <- count + 1
        approxs[[count]] <- approx(x = ind.scan[,i], y = ind.scan[,j], xout = ind.scan[,i])
        xy[i,j] <- count
      } 
    }

  } 
  else if (typ == 2) 
  { # impute from scan numbers
    # obj2 must be a clusterAlignment having not used peaks
	if(is.null(obj2))
	  stop("'obj2' must be a 'clusterAlignment' object (run with usePeaks = F).")
	  	  
    for(i in 1:n){ 
      for(j in 1:n){
        if(i == j) next
	    al <- obj2@alignments[[obj2@aligned[runs[i],runs[j]]]]
        count <- count+1
		if(i < j)
    {
		  xx <- al@v$match[,1]
		  yy <- al@v$match[,2]
		} 
    else
    {
		  xx <- al@v$match[,2]
		  yy <- al@v$match[,1]
		}
        approxs[[count]] <- approx(x = xx, y = yy, xout = ind.scan[,i])
        xy[i,j] <- count
      } 
    }
  }
  
  v <- which(is.na(ind.scan), arr.ind = TRUE)
  v <- v[v[,1] %in% rws,]

  for(i in 1:nrow(v)){
    pr <- xy[,v[i,2]]
    pr <- pr[!is.na(pr)]
    implist <- lapply(approxs[pr], FUN = function(u,ind) u$y[ind], ind = v[i,1])  # use all other runs to impute
    newind[v[i,1],v[i,2]] <- round(median(unlist(implist), na.rm = TRUE))
  }
  
  newind.start <- matrix( ,nrow = nrow(ind), ncol = ncol(ind))
  newind.end <- matrix( ,nrow = nrow(ind), ncol = ncol(ind))
  before <- matrix( ,nrow = nrow(ind), ncol = ncol(ind))
  after <- matrix( ,nrow = nrow(ind), ncol = ncol(ind))
  for(j in 1:ncol(newind.start)){
    before[,j] <- (pD@peaksind[[ runs[j] ]] - pD@peaksind.start[[runs[j]]])[ind[,j]]
    after[,j] <- (pD@peaksind.end[[runs[j]]] - pD@peaksind[[runs[j]]])[ind[,j]]
  }
  before <- round(rowMeans(before, na.rm = TRUE))
  after <- round(rowMeans(after, na.rm = TRUE))
  
  for(i in 1:nrow(newind)){
    for(j in 1:ncol(newind)){
      if(!is.na(newind[i,j]))
      {
	    newind.start[i,j] <- newind[i,j] - before[j]
	    newind.end[i,j] <- newind[i,j] + after[j]
      }
  	}  
  }

  newind.start[newind.start <= 0] <- 1
  
  rownames(newind) <- rownames(newind.start) <- rownames(newind.end) <- 1:nrow(newind)
  list(apex = newind, start = newind.start, end = newind.end)
}


calcTimeDiffs <- function(pd, ca.full, verbose = TRUE){
  n <- length(ca.full@alignments)
  estdiffs <- vector("list", n)
  for(k in 1:n){
    v <- which(ca.full@aligned == k, arr.ind = TRUE)
    run1 <- v[2,1]; run2 <- v[2,2]
	if(run1 > run2)
    {
      tmp <- run2; run2 <- run1; run1 <- tmp
    }
    if(verbose)
	  {
      cat("Estimating time differences for alignment", k, "(", run1, "-", run2,")\n")
    }
    pind1 <- pd@peaksind[[run1]]
    pind2 <- pd@peaksind[[run2]]
    scantimes <- c(median(diff(pd@rawrt[[run1]])), median(diff(pd@rawrt[[run2]])))
    m <- ca.full@alignments[[ca.full@aligned[run1,run2]]]@v$match
    est1 <- approx(m[,1],m[,2], xout = pind1)
    est2 <- approx(m[,2],m[,1], xout = pind2)
    estdiff <- matrix(0, length(pind1), length(pind2))
    for(i in 1:nrow(estdiff)){
      for(j in 1:ncol(estdiff)){
	    diffs <- c(est1$y[i] - est2$x[j], est1$x[i] - est2$y[j])
	    ## cat(i,j,diffs,"\n")
	    am <- which.min(abs(diffs))
	    estdiff[i,j] <- diffs[am] * scantimes[am]
	   }
    }
    estdiffs[[k]] <- estdiff
  }
  estdiffs
}

