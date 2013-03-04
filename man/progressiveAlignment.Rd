\name{progressiveAlignment-class}

\alias{progressiveAlignment-class}
\alias{progressiveAlignment-show}
\alias{progressiveAlignment}
\alias{show,progressiveAlignment-method}

\title{Data Structure for progressive alignment of many GCMS samples}

\description{Performs a progressive peak alignment (clustalw style) of multiple GCMS peak lists}

\usage{progressiveAlignment(pD,cA,D=1000,gap=.5,verbose=TRUE,usePeaks=TRUE,df=30,compress=TRUE)}

\arguments{

\item{pD}{a \code{peaksDataset} object}

\item{cA}{a \code{clusterAlignment} object}

\item{D}{retention time penalty}

\item{gap}{gap parameter}

\item{verbose}{logical, whether to print information}

\item{usePeaks}{logical, whether to use peaks (if \code{TRUE}) or the full 2D profile alignment (if \code{FALSE})}

\item{df}{distance from diagonal to calculate similarity}

\item{compress}{logical, whether to store the similarity matrices in sparse form}

}


\details{

The progressive peak alignment we implemented here for multiple GCMS peak lists is analogous to how \code{clustalw} takes a set of pairwise sequence alignments and progressively builds a multiple alignment.  More details can be found in the reference below.

}

\value{

\code{progressiveAlignment} object

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data \emph{PhD dissertation} University of Melbourne.

}


\seealso{
\code{\link{peaksDataset}}, \code{\link{multipleAlignment}}
}

\examples{
require(gcspikelite)

# paths and files
gcmsPath<-paste(find.package("gcspikelite"),"data",sep="/")
cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
eluFiles<-dir(gcmsPath,"ELU",full=TRUE)

# read data, peak detection results
pd<-peaksDataset(cdfFiles[1:2],mz=seq(50,550),rtrange=c(7.5,8.5))
pd<-addAMDISPeaks(pd,eluFiles[1:2])

ca<-clusterAlignment(pd, gap = .5,D=.05,df=30)
pa<-progressiveAlignment(pd, ca, gap = .6, D=.1,df=30)
}

\keyword{classes}


