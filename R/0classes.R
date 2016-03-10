
setClass("peaksDataset",
#  components of a GCMS dataset
representation(rawdata="list",rawrt="list",mz="numeric",files="character",
peaksdata="list",peaksrt="list",peaksind="list",peaksind.start="list",peaksind.end="list")
)

setClassUnion("eitherMatrix",c("matrix", "matrix.csc"))

setClass("peaksAlignment",
#  components of a pairwise alignment
representation(v="list",r="eitherMatrix",sim="numeric",gap="numeric",dist="numeric",D="numeric",compressed="logical")
)

setGeneric("compress", function(object, verbose=TRUE, ...) standardGeneric("compress"))
setGeneric("decompress", function(object, verbose=TRUE, ...) standardGeneric("decompress"))
setGeneric("plotImage", function(object,run=1,rtrange=c(11,13),main=NULL,mzrange=c(50,200),SCALE=log2,...) standardGeneric("plotImage"))

setClass("clusterAlignment",
representation(runs="integer",aligned="matrix",gap="numeric",D="numeric",dist="matrix",alignments="list",merge="matrix")
)

setClass("progressiveAlignment",
#  components of a progressive Alignment (cluster tree used to guide alignment)
representation(merges="list")
)

setClass("betweenAlignment",
representation(mergedPeaksDataset="peaksDataset",ind="matrix",imputeind="list",
               runs="integer",groups="character",cA="clusterAlignment",pA="progressiveAlignment",filtind="list",newind="list")
)

setClass("multipleAlignment",
# alignment of several samples, possibly with multiple groups
representation(clusterAlignmentList="list",progressiveAlignmentList="list",
timeDiff="list",impute="list",betweenAlignment="betweenAlignment")
)


# setClassUnion("align",   c('matrix',   'data.frame'))
# setClass(Class='correlationAlignment',
#          slots=c(Alignment='align', Center='character'))

setClass(Class='correlationAlignment',
         slots=c(Alignment='data.frame', Center='character'))




