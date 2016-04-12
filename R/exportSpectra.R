
##' Write the deconvoluted mass spectra to an external file
##'
##' Write a .msp file of the deconvoluted mass spectra. Usfull to try
##' to identify the unknown spectra using NIST Search.
##' @title exportSpectra
##' @param object an object of class "peaksDataset" where to keep the
##' mass spectra; both abundance (y) than m/z (x)  
##' @param sample character, the sample from were to plot the mass spectra
##' @param spectraID numerical, a vector containing the index of the spectra to be
##' plotted.
##' @param normalize logical, if TRUE normalize the intensity of the
##' mass peak to 100, the most abundant is 100% and the other peaks
##' are scaled consequetially
##' @return a .msp file ready to be read using NIST search
##' @author riccardo.romoli@unifi.it
exportSpectra <- function(object, sample, spectraID, normalize=TRUE){
    spectraID <- as.numeric(spectraID)
    sample <- as.character(sample)
    i <- grep(pattern=sample, object@files)
    sp <- grep(pattern=spectraID, colnames(object@peaksdata[[i]]))
    if(length(sp) > 0){
        y <- object@peaksdata[[i]][,sp]        
        if(normalize){
            y <- sapply(y, function(x){
                i <- which.max(y)
                (100*x)/y[i]}
                        )}
        spec <- paste(object@mz, y)
        msp <- rbind(paste("NAME: Variable", spectraID),
                     paste("COMMENT:", round(object@peaksrt[[i]][sp], digits=2), 'min'),
                     paste("FORMULA:"),
                     paste("MW:"),
                     paste("CAS:"),
                     paste("SYNONYM:"),
                     paste("Num Peaks:", length(spec)),
                     matrix(unlist(spec), nrow=length(unlist(spec)),
                            ncol=1)
                     )
        # fix the presence of '/' due to full=TRUE in the dir() function
        sample <- sub(pattern='/', replacement='-', sample, perl=TRUE)
        write(msp, file=paste(spectraID, '_', sample, '.msp'), sep="\n")
    }else{
        stop(paste('The spectrum is not present in the sample', sample, '\n'))
    }
}
