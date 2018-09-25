
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
    browser()
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


##' Read the mass spectra from an external msp file
##'
##' Read the mass spectra from an external file in msp format. The format is
##' used in NIST search library database.
##' @title importSpec
##' @param file 
##' @return list conaining the mass spctra
##' @author riccardo.romoli@unifi.it
importSpec <- function(file){
    ## read msp lib
    lib <- scan(file, what = "", sep = "\n", quiet = TRUE)
    ## separate each mass spec
    starts <- which(regexpr("[Nn][Aa][Mm][Ee]: ", lib) == 1)
    ends <- c(starts[-1] - 1, length(lib))
    ## loop to extract the mass spec into a list
    list.spec <- lapply(1:length(starts), function(z){
        ## meta data
        comp <- lib[starts[z]:ends[z]]
        numPeaks.idx <- which(regexpr("Num Peaks: ", comp) == 1)
        metaData <- comp[1:numPeaks.idx - 1]
        md <- strsplit(metaData, split = ": ")
        md1 <- sapply(md, "[[", 1)
        md2 <- sapply(md, "[", 2)
        metaData.list <- setNames(as.list(md2), md1)
        ## mass spec
        nlines <- length(comp)
        npeaks <- as.numeric(strsplit(comp[numPeaks.idx], ":")[[1]][2])
        peaks.idx <- (numPeaks.idx + 1):nlines
        pks <- gsub("^ +", "", unlist(strsplit(comp[peaks.idx], ";")))
        pks <- pks[pks != ""]
        if (length(pks) != npeaks)
            stop("Not the right number of peaks in compound", cmpnd$Name)
        pklst <- strsplit(pks, " ")
        pklst <- lapply(pklst, function(x) x[x != ""])
        mz <- as.numeric(sapply(pklst, "[[", 1))
        int <- as.numeric(sapply(pklst, "[[", 2))
        ## 
        finaltab <- matrix(c(mz, int), ncol = 2)
        if (any(table(mz) > 1))
        {
            warning("Duplicate mass in compound ", cmpnd$Name,
                    " (CAS ", cmpnd$CAS, ")... summing up intensities")
            finaltab <- aggregate(finaltab[,2],
                                  by = list(finaltab[,1]),
                                  FUN = sum)
        }

        colnames(finaltab) <- c("mz", "intensity")
        c(metaData.list, list(spec = finaltab))
        
    }
    )
    return(list.spec)
}
