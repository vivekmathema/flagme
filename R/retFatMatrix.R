#' Build a fat data matrix
#'
#' This function allows to extract the data from an object created using 
#' \code{gatherInfo} and build a data matrix using the area of the deconvoluted
#' and aligned peaks. The row are the samples while the column represent the 
#' different peaks.
#' @title retFatMatrix
#' @param object peakDataset object
#' @param data a gatherInfo() object
#' @param minFilter the minimum number for a feature to be returned 
#' in the data matrix
#'
#' @return A fat data matrix containing the area of the deconvoluted and aligned
#' peaks. The row are the samples while the column represent the different peaks
#' @author Riccardo Romoli \email{riccardo.romoli@unifi.it}
#' @seealso \code{\link{gatherInfo}}
#' @keywords 
#' @examples 
#' require(gcspikelite)
#' # paths and files
#' gcmsPath <- paste(find.package("gcspikelite"), "data", sep = "/")
#' cdfFiles <- dir(gcmsPath,"CDF",full=TRUE)
#' # read data, peak detection results
#' pd <- peaksDataset(cdfFiles[1:2], mz=seq(50,550),
#'                    rtrange=c(7.5,8.5))
#' pd <- addXCMSPeaks(files=cdfFiles[1:2], object=pd,
#'                    peakPicking=c('mF'), snthresh=3, fwhm=4,
#'                    step=1, steps=2, mzdiff=0.5)
#' ma <- multipleAlignment(pd = pd, group = c(1,1),
#'                         filterMin = 1, metric = 2, type = 2)
#' outList <- gatherInfo(pd, ma)
#' mtxD <- retFatMatrix(object = pd, data = outList, minFilter = 1)
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
