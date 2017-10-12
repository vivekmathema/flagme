
#' Plot the mass spectra from the profile matrix 
#'
#' Plot the deconvoluted mass spectra from the profile matrix 
#' @title plotSpectra
#' @param object an object of class "peaksDataset" where to keep the
#' mass spectra; both abundance (y) than m/z (x)  
#' @param sample character, the sample from were to plot the mass spectra
#' @param spectraID numerical, a vector containing the index of the spectra to be
#' plotted.
#' @param normalize logical, if TRUE normalize the intensity of the
#' mass peak to 100, the most abundant is 100% and the other peaks
#' are scaled consequetially
#' @param ... other parameter passed to the plot() function
#' @author riccardo.romoli@unifi.it
#' @examples
#' gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
#' cdfFiles <- dir(gcmsPath,"CDF", full=TRUE)
#' # read data, peak detection results
#' pd <- peaksDataset(cdfFiles[1:3], mz=seq(50,550), rtrange=c(7.5,10.5))
#' pd <- addXCMSPeaks(files=cdfFiles[1:3], object=pd, peakPicking=c('mF'),
#'                    snthresh=3, fwhm=10,  step=0.1, steps=2, mzdiff=0.5,
#'                    sleep=0)
#' ## align two chromatogram
#' pA <- peaksAlignment(pd@peaksdata[[1]], pd@peaksdata[[2]],
#'                      pd@peaksrt[[1]], pd@peaksrt[[2]], D=50,
#'                      metric=3, compress=FALSE, type=2, penality=0.2)
#' pA@v$match
#' ## plot the mass spectra
#' par(mfrow=c(2,1))
#' plotSpectra(object=pd, sample=cdfFiles[1], spectraID=10)
#' plotSpectra(object=pd, sample=cdfFiles[2], spectraID=12)
plotSpectra <- function(object, sample, spectraID, normalize=TRUE, ...){
    ## spectraID <- as.numeric(spectraID)
    sp <- as.numeric(spectraID)
    sample <- as.character(sample)
    i <- grep(pattern=sample, object@files)
    ## sp <- grep(pattern=spectraID, colnames(object@peaksdata[[i]]))
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


#' Plot the aligned mass spectra 
#'
#' Plot the deconvoluted and aligned mass spectra collected using gatherInfo()
#' @title plotMultipleSpectra
#' @param object where to keep the mass range of the experiment
#' @param outList where to keep the mass spectra; both abundance than m/z
#' @param spectra a vector containing the index of the spectra to be plotted.
#'                Is referred to outList
#' @param fullRange if TRUE uses the mass range of the whole experiment, 
#'                  otherwise uses only the mass range of each plotted spectum
#' @param normalize if TRUE normalize the intensity of the mass peak to 100, 
#'                  the most abundant is 100% and the other peaks are scaled 
#'                  consequetially
#' @param ... further arguments passed to the ‘plot’ command
#'
#' @return
#' @author Riccardo Romoli (riccardo.romoli@unifi.it)
#' @seealso
#' @keywords plot() gatherInfo()
#' @examples
#' ## Rd workflow
#' gcmsPath <- paste(find.package("gcspikelite"), "data", sep = "/")
#' cdfFiles <- dir(gcmsPath,"CDF", full = TRUE)
#' 
#' # read data, peak detection results
#' pd <- peaksDataset(cdfFiles[1:4], mz = seq(50,550), rtrange = c(7.5,10.5))
#' pd <- addXCMSPeaks(files = cdfFiles[1:4], object = pd, peakPicking = c('mF'),
#'                    snthresh = 2, fwhm = 8,  step = 0.5, steps = 2, mzdiff = 0.5,
#'                    sleep = 0)
#' ## multiple alignment
#' ma <- multipleAlignment(pd, c(1,1,2,2), wn.gap = 0.5, wn.D = 0.05, bw.gap = 0.6, 
#'                         bw.D = 0.2, usePeaks = TRUE, filterMin = 1, df = 50,
#'                         verbose = TRUE, metric = 2, type = 2)
#' 
#' ## gather apex intensities
#' gip <- gatherInfo(pd, ma)
#' gip[[33]]
#' plotMultipleSpectra(object = pd, outList = gip, spectra = 33, fullRange = FALSE,
#'                     normalize = TRUE)
plotMultipleSpectra <- function(object, outList, spectra, fullRange = TRUE,
                                normalize = TRUE, ...){
spectra <- as.numeric(spectra)
mz <- outList[[spectra]]$mz
abundance <- outList[[spectra]]$data

if(normalize)
{
    for(i in 1:ncol(abundance))
    {
        if(is.na(sum(abundance[,i])))
            next
        abundance[,i] <- 100 * abundance[,i] / abundance[which.max(abundance[,i]), i]
    }
}   

specnum <- table(apply(abundance, 2, sum)) # count number of massSpec different from NA

par(mfrow = c(length(specnum), 1))
for(i in 1:ncol(abundance))
    {
        if(is.na(apply(abundance, 2, sum)[i]))
            next
        plot(mz, abundance[,i], type = 'h', 
            main = paste('Retention Time', 
                         round(mean(outList[[spectra]]$rt[i], na.rm = TRUE), 
                         digits = 2), 'min - Sample ', object@files[i]),
             xlab = 'm/z', 
             ylab = if(normalize)
                    {
                        'Rel. Abundance'
                    } 
                    else 
                    {
                        'Abs. Abundance'
                    }, 
             
             xlim = if(fullRange) 
                    {
                        range(object@mz)
                    }
             )    
    }

}