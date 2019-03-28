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
#' in the data matrix. Default is 2/3 of the samples
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
retFatMatrix <- function (object, data, minFilter = round(length(files)/3*2)) 
{
    a <- lapply(seq(along = data), function(x)
    {
        apply(data[[x]]$data, 2, sum)
    }
    )
    ## i nomi delle colonne equivalgono al numero del file;
    ## il numero della riga equivale alla tasca della lista di gatherInfo()
    abumtx <- do.call(rbind, a) 
    abumtx <- apply(abumtx, 1, "[")
    files_to_merge <- rownames(abumtx)
    files.idx <- as.numeric(sub(pattern = "^.", replacement = "", files_to_merge))
    sample <- object@files[files.idx]
    colnames(abumtx) <- sapply(1:ncol(abumtx), function(x){paste0("Feat", x)})
    mf <- minFilter
    keep <- c()
    for (g in 1:ncol(abumtx))
    {
        keep[g] <- sum(!is.na(abumtx[, g])) >= mf
    }
    abumtx[is.na(abumtx)] <- c(0)
    df <- cbind.data.frame(sample, abumtx[, keep])
    return(df)
}

