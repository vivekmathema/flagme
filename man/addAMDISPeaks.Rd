\name{addAMDISPeaks}
\alias{addAMDISPeaks}

\title{Add AMDIS peak detection results}

\description{Reads ASCII ELU-format files (output from AMDIS) and attaches them to an already created \code{peaksDataset} object}

\usage{addAMDISPeaks(object,fns=dir(,"[Eu][Ll][Uu]"),verbose=TRUE,...)}

\arguments{

\item{object}{a \code{peaksDataset} object.}

\item{fns}{character vector of same length as \code{object@rawdata} (user ensures the order matches)}

\item{verbose}{whether to give verbose output, default \code{TRUE}}

\item{...}{arguments passed on to \code{parseELU}}

}


\details{

Repeated calls to \code{parseELU} to add peak detection results to the original \code{peaksDataset} object.

}

\value{

\code{peaksDataset} object

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data 
\emph{PhD dissertation} University of Melbourne.

}


\seealso{
\code{\link{parseELU}}, \code{\link{peaksDataset}}
}

\examples{
# need access to CDF (raw data) and ELU files 
require(gcspikelite)
gcmsPath<-paste(find.package("gcspikelite"),"data",sep="/")

# full paths to file names
cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
eluFiles<-dir(gcmsPath,"ELU",full=TRUE)

# create a 'peaksDataset' object and add AMDIS peaks to it
pd<-peaksDataset(cdfFiles[1],mz=seq(50,550),rtrange=c(7.5,8.5))
pd<-addAMDISPeaks(pd,eluFiles[1])
}

\keyword{manip}


