retFatMatrix <- function(object, data, minFilter=1){
    a <- lapply(seq(along=data), function(x){
        apply(data[[x]]$data, 2, sum)
    })
    abumtx <- do.call(rbind, a)
    abumtx <- apply(abumtx, 1, '[' ) # works as t()
    ## s <- rownames(abumtx)
    rownames(abumtx) <- c(1:(dim(abumtx)[1]))
    sample <- object@files # sample names
    ## if(length(names(data[[1]]$rt)) == 0)
    ## {
    ##     sample <- data[[1]]$traceRaw[,1] ## ricca
    ## }
    ## else
    ## {
    ##     sample <- names(data[[1]]$rt) ## mark
    ## }
    
    ## min filter
    mf <- minFilter
    keep <- c()
    for(g in 1:ncol(abumtx)){
        keep[g] <- sum(!is.na(abumtx[,g])) >= mf
    }

    abumtx[is.na(abumtx)] <- c(0) # remove NAs
    df <- cbind.data.frame(sample, abumtx[,keep])
    
    return(df)
}
