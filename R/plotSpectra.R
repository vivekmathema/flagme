
##' Plot the mass spectra from the profile matrix 
##'
##' Plot the deconvoluted mass spectra from the profile matrix 
##' @title plotSpectra
##' @param object an object of class "peaksDataset" where to keep the
##' mass spectra; both abundance (y) than m/z (x)  
##' @param sample character, the sample from were to plot the mass spectra
##' @param spectraID numerical, a vector containing the index of the spectra to be
##' plotted.
##' @param normalize logical, if TRUE normalize the intensity of the
##' mass peak to 100, the most abundant is 100% and the other peaks
##' are scaled consequetially
##' @param ... other parameter passed to the plot() function
##' @author riccardo.romoli@unifi.it
##' @examples
##' ## need access to CDF (raw data)
##' require(gcspikelite)
##' gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
##'
##' ## full paths to file names
##' cdfFiles <- dir(gcmsPath, "CDF", full=TRUE)
##'
##' ## create a 'peaksDataset' object and add XCMS peaks to it
##' pd <- peaksDataset(cdfFiles[1:3], mz=seq(50,550), rtrange=c(7.5,8.5))
##' pd <- addXCMSPeaks(cdfFiles[1:3], pd, peakPicking=c('mF'), snthresh=5,
##'                    fwhm=10, step=1, steps=2, mzdiff=0.5, sleep=0)
##' 
##' ## align the chromatograms
##' mp <- correlationAlignment(object=pd, thr=0.8, D=20,
##'                            penality=0.2, normalize=TRUE,
##'                            minFilter=2)
##' ## view the alignment results
##' mp@Alignment
##' 
##' ## plot the mass spectra
##' par(mfrow=c(3,1))
##' plotSpectra(object=pd, sample=cdfFiles[1], spectraID=2)
##' plotSpectra(object=pd, sample=cdfFiles[2], spectraID=3)
##' plotSpectra(object=pd, sample=cdfFiles[3], spectraID=4)
plotSpectra <- function(object, sample, spectraID, normalize=TRUE, ...){
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
        plot(y,
             type='h',
             main=paste('Sample', names(object@peaksdata[i]), 'RT', object@peaksrt[[i]][sp], 'min'),
             xlab='m/z',
             ylab=
                 if(normalize){
                     'Rel. Abundance'
                 }else{
                     'Abs. Abundance'
                 }, ...)
    }else{
        stop(paste('The spectrum is not present in the sample', sample, '\n'))
    }
}
