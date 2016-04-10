retFatMatrix <- function(data){
    a <- lapply(seq(along=data), function(x){
        apply(data[[x]]$data, 2, sum)
    })
    abumtx <- do.call(rbind, a)
    abumtx <- apply(abumtx, 1, '[' ) # works as t()
    s <- rownames(abumtx)
    rownames(abumtx) <- c(1:(dim(abumtx)[1]))
    abumtx[is.na(abumtx)] <- c(0)
    if(length(names(data[[1]]$rt)) == 0){
        sample <- data[[1]]$traceRaw[,1] ## ricca
    }else{
        sample <- names(data[[1]]$rt) ## mark
    }
    df <- cbind.data.frame(sample, abumtx)
    return(df)
}
