
#' A class description
#'
#' @import methods
#' @importClassesFrom SparseM matrix.csc
#' @importMethodsFrom SparseM  as.matrix.csc
#' @importFrom SparseM  as.matrix.csc 
setClassUnion("eitherMatrix", c("matrix", "matrix.csc"))


#' A class description
#'
#' @import methods
#' @export peaksDataset
#' @exportClass peaksDataset
#' @noRd
setClass("peaksDataset",
         ## components of a GCMS dataset
         representation(rawdata="list", rawrt="list", mz="numeric",
                        files="character",
                        peaksdata="list", peaksrt="list",
                        peaksind="list", peaksind.start="list",
                        peaksind.end="list")
         )

setGeneric("plotChrom",
           function(object, runs = 1:length(object@rawdata),
                    mzind = 1:nrow(object@rawdata[[1]]),
                    mind = NULL, plotSampleLabels = TRUE,
                    calcGlobalMax = FALSE, peakCex = 0.8,
                    plotPeaks = TRUE, plotPeakBoundaries = FALSE,
                    plotPeakLabels = FALSE,
                    plotMergedPeakLabels = TRUE, mlwd = 3,
                    usePeaks = TRUE, plotAcrossRuns = FALSE,
                    overlap = F, rtrange = NULL, cols = NULL, thin = 1,
                    max.near = median(object@rawrt[[1]]), how.near = 50,
                    scale.up = 1, ...)
               standardGeneric("plotChrom")
           )


setClass("peaksAlignment",
         ## components of a pairwise alignment
         representation(v="list", r="eitherMatrix", sim="numeric",
                        gap="numeric", dist="numeric", D="numeric",
                        compressed="logical") 
         )

setGeneric("plotAlignment",
           function(object, xlab = "Peaks - run 1",
                    ylab = "Peaks - run 2", plotMatches = TRUE,
                    matchPch = 19, matchLwd = 3, matchCex = 0.5,
                    matchCol = "black",
                    col = colorpanel(50,'white', "green","navyblue"),
                    breaks = seq(0, 1, length = 51), ...)
               standardGeneric("plotAlignment")
           )  



setGeneric("compress", function(object, verbose=TRUE, ...) standardGeneric("compress"))
setGeneric("decompress", function(object, verbose=TRUE, ...) standardGeneric("decompress"))
setGeneric("plotImage",
           function(object, run=1, rtrange=c(11,13), main=NULL,
                    mzrange=c(50,200), SCALE=log2, ...)
    standardGeneric("plotImage"))  



setClass("clusterAlignment",
         representation(runs="integer", aligned="matrix",
         gap="numeric", D="numeric", dist="matrix", alignments="list",
         merge="matrix")
         )

setGeneric("plotClustAlignment",
           function(object, alignment = 1, ...)
               standardGeneric("plotClustAlignment"))  



setClass("progressiveAlignment",
         ## components of a progressive Alignment (cluster tree used to guide alignment)
         representation(merges="list")
         )

setClass("betweenAlignment",
         representation(mergedPeaksDataset="peaksDataset",
         ind="matrix", imputeind="list", runs="integer",
         groups="character", cA="clusterAlignment",
         pA="progressiveAlignment", filtind="list", newind="list")
         )


setClass("multipleAlignment",
## alignment of several samples, possibly with multiple groups
representation(clusterAlignmentList="list",
               progressiveAlignmentList="list", timeDiff="list",
               impute="list", betweenAlignment="betweenAlignment")
)


## setClass(Class='correlationAlignment',
##          slots=c(Alignment='data.frame', Center='character'))
