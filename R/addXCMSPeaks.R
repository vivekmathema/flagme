
addXCMSPeaks <- function(files, object, peakPicking=c('cwt','mF'), ...){
    cdfFiles <- as.character(files)
    if(length(cdfFiles) != length(object@rawdata))
        stop('Number of files must be the same as the number of runs (and must match).')
    ## features extraction and deconvolution
    xs <- lapply(cdfFiles,
                 function(x, y){
                     ## retention time range to scanrange
                     f <- which(cdfFiles %in% x) 
                     xr <- xcmsRaw(x)
                     rtrange <- c(min(object@rawrt[[f]]), max(object@rawrt[[f]]))*60
                     scanRange <- c(max(1,which(xr@scantime > rtrange[1])[1], na.rm=TRUE),
                                    min(length(xr@scantime), which(xr@scantime > rtrange[2])[1] - 1, na.rm=TRUE))
                     ## peak picking
                     if(peakPicking == 'cwt'){
                         s <- xcmsSet(x, method='centWave', peakwidth=c(5,15),
                                      prefilter=c(3,100), scanrange=scanRange, integrate=1,
                                      mzdiff=-0.001, fitgauss=TRUE, ...)
                     }
                     if(peakPicking == 'mF'){
                         s <- xcmsSet(x, method='matchedFilter', scanrange=scanRange,
                                      ## step=0.5, steps=2, mzdiff=0.5,
                                      max=500, ...)
                     }
                     ## deconvolution
                     a <- annotate(s, perfwhm=0.6, max_peaks=500, quick=TRUE) 
                     return(a)
                 }, 
                 y=peakPicking
                 )
    ## which colum to be passed based on the peak-picking step    
    if(peakPicking == 'cwt'){
        area <- c('intb')
    }
    if(peakPicking == 'mF'){
        area <- c('intf')
    }    
    ## build @peaksdata        
    data <- lapply(seq(along = cdfFiles),
                   function(x){
                       filt <- sapply(xs[[x]]@pspectra, function(r){length(r)})
                       spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 6)]
                       mzrange <- object@mz
                       abu <- data.frame(matrix(0, nrow=length(mzrange),
                                                ncol=length(spec.idx)))
                       rownames(abu) <- mzrange
                       colnames(abu) <- spec.idx
                       ## schismatrix
                       mz <- data.frame(mz=mzrange) ## dummy to be merged
                       abu <- sapply(spec.idx, function(z){
                           spec <- getpspectra(xs[[x]], z)[,c('mz', area)]
                           spec[,'mz'] <- round(spec[,'mz'])
                           ## check double mass
                           if (max(table(spec[,1])) > 1){
                               spec.noDouble <- cbind(aggregate(spec[,2],
                                                                list(spec[,1]),
                                                                FUN=sum))
                               colnames(spec.noDouble) <- c('mz', area)
                               spec <- spec.noDouble
                           } else {
                               spec
                           } 
                           abu$z <- merge(spec, mz, by='mz',
                                          all=TRUE)[,area] ## THE GOAL
                       }
                                     )
                       colnames(abu) <- spec.idx ## could be a problem?
                       abu[is.na(abu)] <- c(0)
                       return(abu)
                   }
                   )
    ## build @peaksrt
    apex.rt <- lapply(seq(along=cdfFiles),
                      function(x){
                          filt <- sapply(xs[[x]]@pspectra, function(r){length(r)})
                          spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 6)]
                          apex.rt <- sapply(spec.idx,
                                            function(z){
                                                spec.rt <- getpspectra(xs[[x]], z)[,c('rt')]
                                                rt <- round(mean(spec.rt)/60, digits=3) # in minutes
                                            }
                                            )
                          return(apex.rt)
                      }
                      )    
    ## build @peaksind
    spectra.ind <- lapply(seq(along=cdfFiles),
                          function(x){
                              filt <- sapply(xs[[x]]@pspectra, function(r){length(r)})
                              spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 6)]
                          }
                          )        
    ## build @peaksind.start
    ind.start <- lapply(seq(along=cdfFiles),
                        function(x){
                            filt <- sapply(xs[[x]]@pspectra, function(r){length(r)})
                            spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 6)]
                            rt.start <- sapply(spec.idx,
                                               function(z){
                                                   spec.rt <- getpspectra(xs[[x]], z)[,c('rtmin')]
                                                   rt <- round(mean(spec.rt), digits=3)
                                               }
                                               )
                            return(rt.start)
                        }
                        )
    ## build @peaksind.stop
    ind.stop <- lapply(seq(along=cdfFiles),
                       function(x){
                           filt <- sapply(xs[[x]]@pspectra, function(r){length(r)})
                           spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 6)]
                           rt.stop <- sapply(spec.idx,
                                             function(z){
                                                 spec.rt <- getpspectra(xs[[x]], z)[,c('rtmax')]
                                                 rt <- round(mean(spec.rt), digits=3)
                                             }
                                             )
                           return(rt.stop)
                       }
                       )
    ## build @files
    object@files    
    ## build @mz
    object@mz
    #
    ## S4 ##
    nm <- lapply(files, function(u){sp <- strsplit(u, split="/")[[1]]
                                    sp[length(sp)]
                                }
                 )
    nm <- sub(".CDF$", "" , nm)
    names(data) <- names(apex.rt) <- names(spectra.ind) <- names(ind.start) <- names(ind.stop) <- nm
    
    new("peaksDataset",
        files=object@files,
        peaksdata=data,
        peaksrt=apex.rt,
        peaksind=spectra.ind,
        peaksind.start=ind.start,
        peaksind.end=ind.stop, 
        rawdata=object@rawdata,
        rawrt=object@rawrt,
        mz=object@mz
        )
}

################################################################################
