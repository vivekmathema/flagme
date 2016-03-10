
setMethod("show",
          signature = "correlationAlignment",
          definition = function(object) {
              cat("An object of class ", class(object), "\n", sep='')
              cat('Aligned Samples:', dim(object@Alignment)[2], '\n',
                  'Aligned Peaks:', dim(object@Alignment)[1], '\n',
                  sep='')
              invisible(NULL)
          }
          )


.corP <- function(object, idx1, idx2, D, penality, normalize){
    #
    D <- as.numeric(D) # time window in second
    pn <- as.numeric(penality)# penality if out of time window
    #
    pearson <- function(x,y){
        size <- length(x)
        cfun <- .C("pearson", size=as.integer(size), x=as.double(x),
                   y=as.double(y), result=double(1), PACKAGE='flagme')
        return(cfun[["result"]])
    }
    #
    Normalize <- function(j){
        n <- apply(j, 2, function(k){
            m <- k[which.max(k)]
            norm <- k/m*100
        })
    }
    #
    Rank <- function(u) {
        if (length(u) == 0L) 
            u
        else if (is.matrix(u)) {
            if (nrow(u) > 1L) 
                apply(u, 2L, rank, na.last="keep")
            else row(u)
        }
        else rank(u, na.last="keep")
    }
    #
    if(normalize == TRUE){
        x <- Normalize(object@peaksdata[[idx1]])
        y <- Normalize(object@peaksdata[[idx2]])
    }else if(normalize == FALSE){
        x <- object@peaksdata[[idx1]]
        y <- object@peaksdata[[idx2]]
    }
    method <- c("pearson", "kendall", "spearman")
    ncx <- ncol(x)
    ncy <- ncol(y)
    r <- matrix(0, nrow=ncx, ncol=ncy)
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncy)) {
            x2 <- x[, i]
            y2 <- y[, j]
            ok <- complete.cases(x2, y2)
            x2 <- rank(x2[ok])
            y2 <- rank(y2[ok])
            ## insert rt penality in seconds
            rtDiff <- object@peaksrt[[idx1]][i]*60 - object@peaksrt[[idx2]][j]*60 # retention time in seconds
            rtDiff <- abs(rtDiff)
            r[i, j] <- if (any(ok))
                           if(rtDiff <= D)
                               pearson(x2, y2)
                           else 
                               pearson(x2, y2) - pn
                       else NA
        }
    }
    rownames(r) <- colnames(x)
    colnames(r) <- colnames(y)
    return(r)
}


correlationAlignment <- function(object, thr=0.85, D=20, penality=0.2,
                                 normalize=TRUE, minFilter=1){
    D <- as.numeric(D)
    penality <- as.numeric(penality)
    minFilter <- as.numeric(minFilter)
    #
    if((minFilter > length(object@files)) == TRUE){
        stop(paste('The minimum number for a compound to be matched must not exceed the number of the sample!!!',
                   '\n', 'Please lower the minFil parameter.'))
    }    
    # center sample for alignment
    maxFeat <- max(sapply(object@peaksind, function(x) {length(x)}))# numero massimo picchi 
    center <- names(unlist(sapply(object@peaksind,
                                  function(x){which(length(x) ==
                                                        maxFeat)})))[1]# capita che siano piu di uno
    #
    # calculate correlation matrix; center riga della matrice
    fl <- names(object@peaksind[which(names(object@peaksind) != center)])# elimino center dai giochi
    list.conf <- lapply(fl, function(x){
        ## conf <- cor(object@peaksdata[[center]], object@peaksdata[[x]])
        conf <- .corP(object, idx1=center, idx2=x, D=D,
                     penality=penality, normalize=normalize)
        ww <- which(conf > thr, arr.ind=TRUE) # row and col number
        ww[,1] <- as.numeric(rownames(conf)) # return the spectrum id
        ww[,2] <- as.numeric(colnames(conf)) # return the spectrum id
        colnames(ww) <- c(center, x)
        return(list(conf, ww, x)) # poi levare da return sia x che conf
    })
    # merge the correlation alignment results
    align <- Reduce(function(x, y) merge(x, y, all=TRUE), sapply(list.conf, '[', 2))
    align <- as.data.frame(align)
    #
    # min filter
    mf <- minFilter
    alignNoNA <- align
    alignNoNA[is.na(alignNoNA)] <- c(0)
    zero <- sapply(1:nrow(alignNoNA), function(x){
        sum(alignNoNA[x,] == 0)
    })
    compPres <- sapply(seq(along=zero), function(x){
        dim(alignNoNA)[2] - zero[x]
    })
    align <- align[which(compPres >= mf),]
    
    ## S4 ##
    new('correlationAlignment',
        Alignment=align,
        Center=center)
}
