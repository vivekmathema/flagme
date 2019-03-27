addXCMSPeaks <- function (files, object, peakPicking = c("cwt", "mF"), perfwhm = 0.75, quick = TRUE, ...)
{
    options(warn = -1)
    cdfFiles <- as.character(files)
    if (length(cdfFiles) != length(object@rawdata)) 
        stop("Number of files must be the same as the number of runs (and must match).")
    xs <- lapply(cdfFiles, function(x, y) {
        f <- which(cdfFiles %in% x)
        xr <- xcmsRaw(x)
        rtrange <- c(min(object@rawrt[[f]]), max(object@rawrt[[f]])) * 
            60
        scanRange <- c(max(1, which(xr@scantime > rtrange[1])[1], 
            na.rm = TRUE), min(length(xr@scantime), which(xr@scantime > 
            rtrange[2])[1] - 1, na.rm = TRUE))
        if (peakPicking == "cwt") {
            s <- xcmsSet(x, method = "centWave", prefilter = c(5, 
                100), scanrange = scanRange, integrate = 1, mzdiff = -0.001, 
                fitgauss = TRUE, ...)
        }
        if (peakPicking == "mF") {
            s <- xcmsSet(x, method = "matchedFilter", scanrange = scanRange, 
                max = 500, ...)
        }
        idx <- which(s@peaks[, "mz"] > min(object@mz) & s@peaks[, 
            "mz"] < max(object@mz))
        s@peaks <- s@peaks[idx, ]
        if(quick == TRUE)
        {
            a <- annotate(s, perfwhm = perfwhm, quick = quick)
        }
        if(quick == FALSE)
        {
            a <- annotate(s, perfwhm = perfwhm, cor_eic_th = 0.8,
                          pval = 0.05, graphMethod = "hcs",
                          calcIso = FALSE, calcCiS = TRUE,
                          calcCaS = FALSE)
        }
        return(a)
    }, y = peakPicking)
    if (peakPicking == "cwt") {
        area <- c("intb")
    }
    if (peakPicking == "mF") {
        area <- c("intf")
    }
    data <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 
            6)]
        mzrange <- object@mz
        abu <- data.frame(matrix(0, nrow = length(mzrange), ncol = length(spec.idx)))
        rownames(abu) <- mzrange
        colnames(abu) <- spec.idx
        mz <- data.frame(mz = mzrange)
        abu <- sapply(spec.idx, function(z) {
            spec <- getpspectra(xs[[x]], z)[, c("mz", area)]
            spec[, "mz"] <- round(spec[, "mz"])
            if (max(table(spec[, 1])) > 1) {
                spec.noDouble <- cbind(aggregate(spec[, 2], list(spec[, 
                  1]), FUN = sum))
                colnames(spec.noDouble) <- c("mz", area)
                spec <- spec.noDouble
            }
            else {
                spec
            }
            abu$z <- merge(spec, mz, by = "mz", all = TRUE)[, 
                area]
        })
        colnames(abu) <- spec.idx
        abu[is.na(abu)] <- c(0)
        return(abu)
    })
    apex.rt <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 
            6)]
        apex.rt <- sapply(spec.idx, function(z) {
            spec.rt <- getpspectra(xs[[x]], z)[, c("rt")]
            rt <- round(mean(spec.rt)/60, digits = 3)
        })
        return(apex.rt)
    })
    spectra.ind <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 
            6)]
    })
    ind.start <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 
            6)]
        rt.start <- sapply(spec.idx, function(z) {
            spec.rt <- getpspectra(xs[[x]], z)[, c("rtmin")]
            rt <- round(mean(spec.rt), digits = 3)
        })
        return(rt.start)
    })
    ind.stop <- lapply(seq(along = cdfFiles), function(x) {
        filt <- sapply(xs[[x]]@pspectra, function(r) {
            length(r)
        })
        spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 
            6)]
        rt.stop <- sapply(spec.idx, function(z) {
            spec.rt <- getpspectra(xs[[x]], z)[, c("rtmax")]
            rt <- round(mean(spec.rt), digits = 3)
        })
        return(rt.stop)
    })
    object@files
    object@mz
    for (i in 1:length(files)) {
        ord <- order(apex.rt[[i]])
        data[[i]] <- data[[i]][, ord]
        apex.rt[[i]] <- apex.rt[[i]][ord]
        spectra.ind[[i]] <- spectra.ind[[i]][ord]
        ind.start[[i]] <- ind.start[[i]][ord]
        ind.stop[[i]] <- ind.stop[[i]][ord]
    }
    options(warn = 0)
    nm <- lapply(files, function(u) {
        sp <- strsplit(u, split = "/")[[1]]
        sp[length(sp)]
    })
    nm <- sub(".CDF$", "", nm)
    names(data) <- names(apex.rt) <- names(spectra.ind) <- names(ind.start) <- names(ind.stop) <- nm
    new("peaksDataset", files = object@files, peaksdata = data, 
        peaksrt = apex.rt, peaksind = spectra.ind, peaksind.start = ind.start, 
        peaksind.end = ind.stop, rawdata = object@rawdata, rawrt = object@rawrt, 
        mz = object@mz)
}



addXCMSPeaks2 <- function(files, object, param = NULL, ...){
    options(warn = -1) # Warning :implicit list embedding of S4 objects is deprecated
    cdfFiles <- as.character(files)
    if(length(cdfFiles) != length(object@rawdata))
        stop('Number of files must be the same as the number of runs (and must match).')
    ## features extraction and deconvolution
    xs <- lapply(cdfFiles,
                 function(x)
                 {
                     ## retention time range to scanrange
                     f <- which(cdfFiles %in% x) 
                     xr <- xcmsRaw(x)
                     ## new xcms
                     raw_data <- readMSData(files = x,  mode = "onDisk")
                                          
                     ## xr@scantime BECAME raw_data@featureData@data$retentionTime
                     rtrange <- c(min(object@rawrt[[f]]),
                                  max(object@rawrt[[f]])) * 60
                     ## scanRange <- c(max(1,
                     ##                    which(xr@scantime > rtrange[1])[1], na.rm = TRUE),
                     ##                min(length(xr@scantime),
                     ##                    which(xr@scantime > rtrange[2])[1] - 1, na.rm = TRUE))
                     scanRange <- c(max(1,
                                        which(raw_data@featureData@data$retentionTime > rtrange[1])[1], na.rm = TRUE),
                                    min(length(raw_data@featureData@data$retentionTime),
                                        which(raw_data@featureData@data$retentionTime > rtrange[2])[1] - 1, na.rm = TRUE))
                     
                     ## peak picking
                     if(class(param) == "NULL")
                     {
                         cat("Please select the integration parametere related to the chromatographic peak detection you want to use \n")
                     }
                     xs <- findChromPeaks(raw_data, param = param)
                     ## We can also coerce the XCMSnExp object into an xcmsSet object:
                     s <- as(xs, "xcmsSet")
                     ## browser()
                     ## rts <- sapply(s@rt$raw, function(x){round(x)})
                     ## mydata <- rts
                    
                     ## require(mclust)
                     ## d_clust <- Mclust(mydata, G = 1:15)
                     ## ## d_clust$BIC
                     ## ## summary(d_clust)
                     ## ##  plot(d_clust)
                     ## rtgroups <- d_clust$classification
                     ## dfrt <- data.frame(feat = names(rtgroups), grp = rtgroups)
                     ## xy.list <- split(dfrt$feat, dfrt$grp)
                     ## xy.list <- lapply(xy.list, function(x){
                     ##     as.numeric(gsub("^F1.S", replacement = "", x = x))})                    
                     ## provo a sostituire lo slot @pspectra con drft
                     
                     ## set mz range
                     idx <- which(
                         s@peaks[,"mz"] > min(object@mz) &
                         s@peaks[,"mz"] < max(object@mz))
                     s@peaks <- s@peaks[idx,]
                     ## deconvolution
                     a <- annotate(s,...)
                     ## a@pspectra <- xy.list
                     ## a <- xsAnnotate(s)
                     ## a <- groupFWHM(a, perfwhm = 0.75)
                     ## a <- groupCorr(a, cor_eic_th = 0.8,
                     ##                pval = 0.05,
                     ##                graphMethod = "hcs",
                     ##                calcIso = FALSE,
                     ##                calcCiS = TRUE,
                     ##                calcCaS = FALSE,
                     ##                cor_exp_th = 0.85) #cor_exp_th = 0.75
                     return(a)
                 }
                 )
    ## which colum to be passed based on the peak-picking step    
    if(class(param)[1] == "CentWaveParam"){
        area <- c('intb')
    }
    if(class(param)[1] == "MatchedFilterParam"){
        area <- c('intf')
    }    
    ## build @peaksdata        
    data <- lapply(seq(along = cdfFiles),
                   function(x){
                       filt <- sapply(xs[[x]]@pspectra, function(r){length(r)})
                       spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 6)]
                       mzrange <- object@mz
                       abu <- data.frame(matrix(0, nrow = length(mzrange),
                                                ncol = length(spec.idx)))
                       rownames(abu) <- mzrange
                       colnames(abu) <- spec.idx
                       ## schismatrix
                       mz <- data.frame(mz = mzrange) ## dummy to be merged
                       abu <- sapply(spec.idx, function(z){
                           spec <- getpspectra(xs[[x]], z)[,c('mz', area)]
                           spec[,'mz'] <- round(spec[,'mz'])
                           ## check double mass
                           if (max(table(spec[,1])) > 1){
                               spec.noDouble <- cbind(aggregate(spec[,2],
                                                                list(spec[,1]),
                                                                FUN = sum))
                               colnames(spec.noDouble) <- c('mz', area)
                               spec <- spec.noDouble
                           } else {
                               spec
                           } 
                           ## bug due to a different mz range from
                           ## peaksDataset() and addXCMS(); SOLVED
                           abu$z <- merge(spec, mz, by = 'mz',
                                          all = TRUE)[,area] ## THE GOAL
                       }
                                     )
                       colnames(abu) <- spec.idx ## could be a problem?
                       abu[is.na(abu)] <- c(0)
                       return(abu)
                   }
                   )
    ## build @peaksrt
    apex.rt <- lapply(seq(along = cdfFiles),
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
    spectra.ind <- lapply(seq(along = cdfFiles),
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
    ind.stop <- lapply(seq(along = cdfFiles),
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
    
    ## order peaks according to retention time
    for(i in 1:length(files)){
        ord <- order(apex.rt[[i]])
        data[[i]] <- data[[i]][,ord]
        apex.rt[[i]] <- apex.rt[[i]][ord]
        spectra.ind[[i]] <- spectra.ind[[i]][ord]
        ind.start[[i]] <- ind.start[[i]][ord]
        ind.stop[[i]] <- ind.stop[[i]][ord]
    }

    options(warn = 0)

    
    ## S4 ##
    nm <- lapply(files, function(u){sp <- strsplit(u, split = "/")[[1]]
                                    sp[length(sp)]
                                }
                 )
    nm <- sub(".CDF$", "" , nm)
    names(data) <- names(apex.rt) <- names(spectra.ind) <- names(ind.start) <- names(ind.stop) <- nm

    new("peaksDataset",
        files = object@files,
        peaksdata = data,
        peaksrt = apex.rt,
        peaksind = spectra.ind,
        peaksind.start = ind.start,
        peaksind.end = ind.stop, 
        rawdata = object@rawdata,
        rawrt = object@rawrt,
        mz = object@mz
        )
}


################################################################################
.flatten3 <- function(x){
    repeat{
        if(!any(vapply(x, is.list, logical(1))))
            return(x)
        x <- Reduce(c, x)
    }
}


.FWHM <- function(a){
    d <- spline(a)
    xmax <- d$x[d$y==max(d$y)]
    x1 <- d$x[d$x < xmax][which.min(abs(d$y[d$x < xmax]-max(d$y)/2))]
    x2 <- d$x[d$x > xmax][which.min(abs(d$y[d$x > xmax]-max(d$y)/2))]
    ## plot(d, type = "l")
    ## points(c(x1, x2), c(d$y[d$x==x1], d$y[d$x==x2]), col="red")
    FWHM <- x2-x1
    return(FWHM)
}


.bayesDeconvolution <- function(xdata, raw_data){
    rts <- unique(xdata[,"rt"])
    dec <- list()
    for(i in seq(rts)){
        idx <-  which(round(xdata[,"rt"], digits = 3) == rts[i])
        if(length(idx) == 0)
        {
            dec[[i]] <- FALSE
        }
        else
        {
            dec[[i]] <- xdata[idx,]
        }
    }
    
    whichRowToTake <- unlist(sapply(dec, function(x){
                length(x) >= 8*12
            }))

    bo <- list()
    for(j in seq(along = dec))
    {
        if(whichRowToTake[j] == FALSE)
        {
            bo[[j]] <- NULL
        }
        else
        {
            bo[[j]] <- dec[[j]]
        }
    }

    # cleanup
    bo[sapply(bo, is.null)] <- NULL

    ## calculate and extract fwhm of the peaks
    ## parallel
    ## clust <- makeCluster(detectCores())
    ## clusterExport(cl = clust, varlist = c(".FWHM", "xdata", "raw_data", "bo"), envir = environment())
    ## fwhm.list <- parLapply(clust, seq(along = bo), function(x){
    fwhm.list <- lapply(seq(along = bo), function(x){
        dec.list <- bo[[x]]
        fwhm.vec <- sapply(1:nrow(dec.list), function(k){
            chr <- chromatogram(raw_data, mz = c(dec.list[k,"mz"], dec.list[k,"mz"]+1),
                                rt = dec.list[,c("rtmin", "rtmax")])
            int <- intensity(chr@.Data[[1]])
            int[is.na(int)] <- 0
            fwhm.ret <- round(.FWHM(int), digits = 3)
            fwhm.ret <- ifelse(test = length(fwhm.ret) == 0, yes = 0, no = fwhm.ret)
            return(fwhm.ret)
        })
        require(mclust)
        d_clust <- mclust::Mclust(fwhm.vec, G = 1:15)
        spec.df <- cbind.data.frame(dec.list, grp = d_clust$classification)
        grp.list <- split(spec.df, spec.df[,"grp"])
        ## clean list
        grp.list.clean <- lapply(grp.list, function(x){
            if(nrow(x) < 6)
            {
                x <- NULL
            }
            else
            {
                x
            }
        })
        ## remove NULL folder
        grp.list.clean[sapply(grp.list.clean, is.null)] <- NULL
        if(length(grp.list.clean) == 0)
        {
            grp.list.clean <- NULL
        }
        else
        {
            names(grp.list.clean) <- paste(x, "-", seq(along = names(grp.list.clean)), sep = "")
        }
        return(grp.list.clean)

    }## end lapply function
    )
    ## stopCluster(clust)
    return(fwhm.list)
}


addXCMSPeaks3 <- function(files, object, param = NULL, ...){
    options(warn = -1) # Warning :implicit list embedding of S4 objects is deprecated
    cdfFiles <- as.character(files)
    if(length(cdfFiles) != length(object@rawdata))
        stop('Number of files must be the same as the number of runs (and must match).')

    ## parallel
    ## clust <- makeCluster(detectCores())
    ## clusterExport(cl = clust, varlist = c("files", "object", "cdfFiles", "param", ".bayesDeconvolution", ".FWHM", ".flatten3"), envir = environment())
    ## xs <- parLapply(clust, cdfFiles, function(x){
    xs <- lapply(cdfFiles, function(x){
    ## features extraction and deconvolution
    ## xs <- lapply(cdfFiles, function(x){
        ## retention time range to scanrange
        f <- which(cdfFiles %in% x) 
        ## xr <- xcmsRaw(x)
        ## new xcms interface
        raw_data <- readMSData(files = x,  mode = "onDisk")        
        ## xr@scantime BECAME raw_data@featureData@data$retentionTime
        rtrange <- c(min(object@rawrt[[f]]),
                     max(object@rawrt[[f]])) * 60
        scanRange <- c(max(1,
                           which(raw_data@featureData@data$retentionTime > rtrange[1])[1], na.rm = TRUE),
                       min(length(raw_data@featureData@data$retentionTime),
                           which(raw_data@featureData@data$retentionTime > rtrange[2])[1] - 1, na.rm = TRUE))                     
        ## peak picking
        if(class(param) == "NULL")
        {
            cat("Please select the integration parametere related to the chromatographic peak detection you want to use \n")
        }
        xs <- findChromPeaks(raw_data, param = param)        
        ## deconvolution
        xd <- chromPeaks(xs)
        a <- .bayesDeconvolution(xdata = xd, raw_data = raw_data)                     
        flat.grp.list <- .flatten3(a)
        decon.list <- sapply(seq(from = 1, to = length(flat.grp.list), by = 13), function(i){
            data.frame(flat.grp.list[i:(i+12)])
        }, simplify = FALSE)
        
        return(decon.list)
    }
    )
    ## stopCluster(clust)
    
    ## which colum to be passed based on the peak-picking step
    if(class(param)[1] == "CentWaveParam")
    {
        area <- c('intb')
    }
    if(class(param)[1] == "MatchedFilterParam")
    {
        area <- c('intf')
    }
    
    ## build @peaksdata
    data <- lapply(seq(along = cdfFiles), function(x){
        ## getpspectra return a matrix like chrompeaks
        ## filt <- sapply(xs[[x]]@pspectra, function(r){length(r)})
        ## spec.idx <- c(1:length(xs[[x]]@pspectra))[which(filt >= 6)]
        spec.idx <- c(1:length(xs[[x]]))
        mzrange <- object@mz
        abu <- data.frame(matrix(0, nrow = length(mzrange),
                                 ncol = length(spec.idx)))
        rownames(abu) <- mzrange
        colnames(abu) <- spec.idx
        ## schismatrix
        mz <- data.frame(mz = mzrange) ## dummy to be merged
        abu <- sapply(spec.idx, function(z){
            spec <- xs[[x]][[z]][,c('mz', area)]# return a matrix
            spec[,'mz'] <- round(spec[,'mz'])
            ## check double mass
            if(max(table(spec[,1])) > 1)
            {
                spec.noDouble <- cbind(aggregate(spec[,2], list(spec[,1]), FUN = sum))
                colnames(spec.noDouble) <- c('mz', area)
                spec <- spec.noDouble
            }
            else
            {
                spec
            } 
            ## bug due to a different mz range from
            ## peaksDataset() and addXCMS(); SOLVED
            abu$z <- merge(spec, mz, by = 'mz', all = TRUE)[,area] ## THE GOAL
        }
        )
        colnames(abu) <- spec.idx ## could be a problem?
        abu[is.na(abu)] <- c(0)
        return(abu)
    }
    )
    
    ## build @peaksrt
    apex.rt <- lapply(seq(along = cdfFiles), function(x){
        spec.idx <- c(1:length(xs[[x]]))
        apex.rt <- sapply(spec.idx, function(z){
            spec.rt <- xs[[x]][[z]][,c("rt")]# return a matrix
            rt <- round(mean(spec.rt)/60, digits = 3) # in minutes
        }
        )
        return(apex.rt)
    }
    )
    
    ## build @peaksind
    spectra.ind <- lapply(seq(along = cdfFiles), function(x){
        spec.idx <- c(1:length(xs[[x]]))
    }
    )

    ## build @peaksind.start
    ind.start <- lapply(seq(along = cdfFiles), function(x){
        spec.idx <- c(1:length(xs[[x]]))
        rt.start <- sapply(spec.idx, function(z){
            spec.rt <- xs[[x]][[z]][,c("rtmin")]
            rt <- round(mean(spec.rt), digits = 3)
        }
        )
        return(rt.start)
    }
    )
    
    ## build @peaksind.stop
    ind.stop <- lapply(seq(along = cdfFiles), function(x){
        spec.idx <- c(1:length(xs[[x]]))
        rt.stop <- sapply(spec.idx, function(z){
            spec.rt <- xs[[x]][[z]][,c("rtmax")]
            rt <- round(mean(spec.rt), digits = 3)
        }
        )
        return(rt.stop)
    }
    )

    ## build @files
    object@files    

    ## build @mz
    object@mz
    
    ## order peaks according to retention time
    for(i in 1:length(files)){
        ord <- order(apex.rt[[i]])
        data[[i]] <- data[[i]][,ord]
        apex.rt[[i]] <- apex.rt[[i]][ord]
        spectra.ind[[i]] <- spectra.ind[[i]][ord]
        ind.start[[i]] <- ind.start[[i]][ord]
        ind.stop[[i]] <- ind.stop[[i]][ord]
    }

    options(warn = 0)

    
    ## S4 ##
    nm <- lapply(files, function(u){sp <- strsplit(u, split = "/")[[1]]
                                    sp[length(sp)]
                                }
                 )
    nm <- sub(".CDF$", "" , nm)
    names(data) <- names(apex.rt) <- names(spectra.ind) <- names(ind.start) <- names(ind.stop) <- nm

    new("peaksDataset",
        files = object@files,
        peaksdata = data,
        peaksrt = apex.rt,
        peaksind = spectra.ind,
        peaksind.start = ind.start,
        peaksind.end = ind.stop, 
        rawdata = object@rawdata,
        rawrt = object@rawrt,
        mz = object@mz
        )
}
