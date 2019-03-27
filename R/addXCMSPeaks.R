addXCMSPeaks <- function (files, object, peakPicking = c("cwt", "mF"), perfwhm = 0.75, quick = quick, ...)
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
