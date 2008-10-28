\name{plotImage}

\alias{plotImage}
\alias{plotImage,peaksDataset-method}

\title{Plot of images of GCMS data}

\description{Image plots (i.e. 2D heatmaps) of raw GCMS profile data}

\usage{plotImage(object,run=1,rtrange=c(11,13),main=NULL,mzrange=c(50,200),SCALE=log2,...)}

\arguments{

\item{object}{a \code{peaksDataset} object}

\item{run}{index of the run to plot an image for}

\item{rtrange}{vector of length 2 giving start and end of the X-axis (retention time)}

\item{main}{main title (auto-constructed if not specified)}

\item{mzrange}{vector of length 2 giving start and end of the Y-axis (mass-to-charge ratio)}

\item{SCALE}{function called to scale the data (default: \code{log2})}

\item{...}{further arguments passed to the \code{image} command}

}


\details{
For \code{peakDataset} objects, each TIC is scale to the maximum value (as specified by the \code{how.near} and \code{max.near} values).  The many parameters gives considerable flexibility of how the TICs can be visualized.

For \code{peakAlignment} objects, the similarity matrix is plotted and optionally, the set of matching peaks.  \code{clusterAlignment} objects are just a collection of all pairwise \code{peakAlignment} objects.
}

\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data \emph{PhD dissertation} University of Melbourne.

}

\seealso{
\code{\link{plot}}, \code{\link{peaksDataset}}
}


\examples{
require(gcspikelite)

# paths and files
gcmsPath<-paste(.find.package("gcspikelite"),"data",sep="/")
cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
eluFiles<-dir(gcmsPath,"ELU",full=TRUE)

# read data
pd<-peaksDataset(cdfFiles[1],mz=seq(50,550),rtrange=c(7.5,8.5))

# image plot
plotImage(pd,run=1,rtrange=c(7.5,8.5),main="")
}

\keyword{classes}


