\name{compress}
\alias{compress}
\alias{decompress}

\alias{compress,clusterAlignment-method}
\alias{decompress,clusterAlignment-method}
\alias{compress,peaksAlignment-method}
\alias{decompress,peaksAlignment-method}
\alias{compress,progressiveAlignment-method}
\alias{decompress,progressiveAlignment-method}


\title{Compress an alignment object}

\description{Many of the peaks are not similar.  So, the set of pairwise similarity matrices can be compressed.}

\usage{
compress(object,verbose=TRUE,...)
decompress(object,verbose=TRUE,...)
}

\arguments{

\item{object}{a \code{peaksAlignment}, \code{peaksAlignment} or \code{peaksAlignment} object to be compressed}

\item{verbose}{logical, whether to print out information}

\item{...}{further arguments}

}


\details{

Using sparse matrix representations, a significant compression can be achieved.  Here, we use the \code{matrix.csc} class of the \code{SpareM} package.

}

\value{

an object of the same type as the input object

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data 
\emph{PhD dissertation} University of Melbourne.

}


\seealso{
\code{\link{peaksAlignment}}, \code{\link{clusterAlignment}}, \code{\link{progressiveAlignment}}
}

\examples{
require(gcspikelite)

# paths and files
gcmsPath<-paste(.find.package("gcspikelite"),"data",sep="/")
cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
eluFiles<-dir(gcmsPath,"ELU",full=TRUE)

# read data, peak detection results
pd<-peaksDataset(cdfFiles[1:2],mz=seq(50,550),rtrange=c(7.5,8.5))
pd<-addAMDISPeaks(pd,eluFiles[1:2])

# pairwise alignment (it is compressed by default)
ca<-clusterAlignment(pd, usePeaks = TRUE, df = 20)
object.size(ca)

# decompress 
ca<-decompress(ca)
object.size(ca)
}

\keyword{manip}


