retFatMatrix <- function(data, minFilter=1){
    a <- lapply(seq(along=data), function(x){
        apply(data[[x]]$data, 2, sum)
    })
    abumtx <- do.call(rbind, a)
    abumtx <- apply(abumtx, 1, '[' ) # works as t()
    s <- rownames(abumtx)
    rownames(abumtx) <- c(1:(dim(abumtx)[1]))
    abumtx[is.na(abumtx)] <- c(0)
    if(length(names(data[[1]]$rt)) == 0)
    {
        sample <- data[[1]]$traceRaw[,1] ## ricca
    }
    else
    {
        sample <- names(data[[1]]$rt) ## mark
    }
    
    ## min filter
    mf <- minFilter
    zerosAllowed <- dim(abumtx)[1]-mf
    ret <- apply(abumtx, 2, table)
    zeros <- lapply(ret, function(x){
        oss <- as.numeric(x[which(as.numeric(names(x)) == 0)])
    }) # count the number of zeros
    ok <- as.numeric(which(zeros <= zerosAllowed)) # the filter
    df <- cbind.data.frame(sample, abumtx[,ok])
    
    return(df)
}
