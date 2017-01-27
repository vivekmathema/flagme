\name{peaksAlignment-class}

\alias{peaksAlignment-class}
\alias{peaksAlignment-show}
\alias{peaksAlignment-plot}
\alias{peaksAlignment}
\alias{show,peaksAlignment-method}
\alias{plot,peaksAlignment-method}
\alias{plot,peaksAlignment,ANY-method}

\title{Data Structure for pairwise alignment of 2 GCMS samples}

\description{Store the raw data and optionally, information regarding
  signal peaks for a number of GCMS runs} 

\usage{peaksAlignment(d1, d2, t1, t2, gap=0.5, D=50, timedf=NULL, df=30,
               verbose=TRUE, usePeaks=TRUE, compress=TRUE, metric=2,
               type=2)}

\arguments{
  \item{d1}{matrix of MS intensities for 1st sample (if doing a peak
    alignment, this contains peak apexes/areas; if doing a profile
    alignment, this contains scan intensities.  Rows are m/z bins, columns
    are peaks/scans.} 
  \item{d2}{matrix of MS intensities for 2nd sample}
  \item{t1}{vector of retention times for 1st sample}
  \item{t2}{vector of retention times for 2nd sample}
  \item{gap}{gap penalty for dynamic programming algorithm. Not used if \code{type=2}}
  \item{D}{time window (on same scale as retention time differences,
    \code{t1} and \code{t2}. Default scale is seconds.)} %% time penalty
  \item{timedf}{list (length = the number of pairwise alignments) of
    matrices giving the expected time differences expected at each pair of
    peaks used with \code{usePeaks}=\code{TRUE}.}
  \item{df}{integer, how far from the diagonal to go to calculate the
    similarity of peaks. Smaller value should run faster, but be careful
    not to choose too low.} 
  \item{verbose}{logical, whether to print out info.}
  \item{usePeaks}{logical, \code{TRUE} uses \code{peakdata} list,
    \code{FALSE} uses \code{rawdata} list for computing similarity.}
  \item{compress}{logical, whether to compress the similarity matrix into
    a sparse format.}
  \item{metric}{numeric, different algorithm to calculate the similarity
  matrix between two mass spectrum. \code{metric=1} call \code{normDotProduct()};
  \code{metric=2} call \code{ndpRT()}; \code{metric=3} call \code{corPrt()}}
  \item{type}{numeric, two different type of alignment function}
}

\details{ 
  peaksAlignment is a hold-all data structure of the raw and peak
  detection data. 
}

\value{
  \code{peaksAlignment} object
}

\author{Mark Robinson, Riccardo Romoli}

\references{
  Mark D Robinson (2008).  Methods for the analysis of gas
  chromatography - mass spectrometry data \emph{PhD dissertation}
  University of Melbourne.  
}


\seealso{
  \code{\link{peaksDataset}}, \code{\link{clusterAlignment}}
}

\examples{
## see clusterAlignment, it calls peaksAlignment

## Not Run:
gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
cdfFiles <- dir(gcmsPath,"CDF", full=TRUE)

# read data, peak detection results
pd <- peaksDataset(cdfFiles[1:3], mz=seq(50,550), rtrange=c(7.5,10.5))
pd <- addXCMSPeaks(files=cdfFiles[1:3], object=pd, peakPicking=c('mF'),
                   snthresh=3, fwhm=10,  step=0.1, steps=2, mzdiff=0.5,
                   sleep=0)
## review peak picking
plot(pd, rtrange=c(7.5, 10.5), runs=c(1:3))

## align two chromatogram
pA <- peaksAlignment(pd@peaksdata[[1]], pd@peaksdata[[2]],
                     pd@peaksrt[[1]], pd@peaksrt[[2]], D=50,
                     metric=3, compress=FALSE, type=2)

plot(pA)
pA@v$match

par(mfrow=c(2,1))
plot(pd@peaksdata[[1]][,15], type='h', main=paste(pd@peaksrt[[1]][[15]]))
plot(pd@peaksdata[[2]][,17], type='h',
     main=paste(pd@peaksrt[[2]][[17]]))
## End (Not Run)
}

\keyword{classes}


