matchSpec <- function(spec1, outList, whichSpec){
    ## first get the average spec from the gatherList
    averageInt <- apply(outList[[whichSpec]]$data, MARGIN = 1, FUN = mean, na.rm = TRUE)
    ## normalize the intensity
    normInt <- sapply(averageInt, function(x){
        i <- which.max(averageInt)
        (100*x)/averageInt[i]
    }
    )
    ## combine mz and normalized intensity
    pspec <- matrix(c(outList[[whichSpec]]$mz, normInt), ncol = 2)
    colnames(pspec) <- c("mz", "intensity")

    libSpec <- spec1$spec

    ## merge the spectra to the same mz range and remove the NA
    mztomerge <- data.frame(mz = min(c(min(libSpec[,"mz"], na.rm = TRUE), min(pspec[,"mz"], na.rm = TRUE))) :
                                max(c(max(libSpec[,"mz"], na.rm = TRUE), max(pspec[,"mz"], na.rm = TRUE)))
                            )
    pspec.merged <- merge(pspec, mztomerge, by = "mz", all = TRUE)
    libSpec.merged <- merge(libSpec, mztomerge, by = "mz", all = TRUE)
    pspec.merged[is.na(pspec.merged)] <- 0
    libSpec.merged[is.na(libSpec.merged)] <- 0
    ## calculate the distance among the spectra
    distance <- normDotProduct(as.matrix(pspec.merged[,"intensity"]), as.matrix(libSpec.merged[,"intensity"]))
    return(distance)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param specFromLib 
##' @param specFromList 
##' @return 
##' @author r
headToTailPlot <- function(specFromLib, specFromList){
    libSpec <- specFromLib$spec
    pspec <- specFromList
    
    ## get average and normalized intensity of the mass spec
    averageInt <- apply(pspec$data, MARGIN = 1, FUN = mean, na.rm = TRUE)
    ## normalize the intensity
    normInt <- sapply(averageInt, function(x){
        i <- which.max(averageInt)
        (1000*x)/averageInt[i]
    }
    )
    pspec.av <- data.frame(mz = pspec$mz, intensity = normInt)

    ## merge the spectra to the same mz range and remove the NA
    mztomerge <- data.frame(mz = min(c(min(libSpec[,"mz"], na.rm = TRUE), min(pspec.av[,"mz"], na.rm = TRUE))) :
                                max(c(max(libSpec[,"mz"], na.rm = TRUE), max(pspec.av[,"mz"], na.rm = TRUE)))
                            )
    pspec.merged <- merge(pspec.av, mztomerge, by = "mz", all = TRUE)
    libSpec.merged <- merge(libSpec, mztomerge, by = "mz", all = TRUE)
    pspec.merged[is.na(pspec.merged)] <- 0
    libSpec.merged[is.na(libSpec.merged)] <- 0

    ## now the plot
    plot(pspec.merged, type = "h", ylim = c(-1000, 1000), main = "Head to Tail Plot")
    points(libSpec.merged[,"mz"], y = -libSpec.merged[,"intensity"], col = 2, type = "h")    
}
