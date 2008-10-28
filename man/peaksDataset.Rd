\name{peaksDataset}

\alias{peaksDataset-class}
\alias{peaksDataset-show}
\alias{peaksDataset-plot}
\alias{peaksDataset}
\alias{show,peaksDataset-method}
\alias{plot,peaksDataset-method}



\title{Data Structure for raw GCMS data and peak detection results}

\description{Store the raw data and optionally, information regarding signal peaks for a number of GCMS runs}

\usage{peaksDataset(fns=dir(,"[Cc][Dd][Ff]"),verbose=TRUE,mz=seq(50,550),rtDivide=60,rtrange=NULL)}

\arguments{

\item{fns}{character vector, filenames of raw data in CDF format.}

\item{verbose}{logical, if \code{TRUE} then iteration progress information is output.}

\item{mz}{vector giving bins of raw data table.}

\item{rtDivide}{number giving the amount to divide the retention times by.}

\item{rtrange}{retention time range to limit data to (must be \code{numeric} vector of length 2)}

}


\details{

peaksDataset is a hold-all data structure of the raw and peak detection data.

}

\value{

\code{peaksDataset} object

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data \emph{PhD dissertation} University of Melbourne.

}


\examples{
require(gcspikelite)

# paths and files
gcmsPath<-paste(.find.package("gcspikelite"),"data",sep="/")
cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
eluFiles<-dir(gcmsPath,"ELU",full=TRUE)

# read data
pd<-peaksDataset(cdfFiles[1:2],mz=seq(50,550),rtrange=c(7.5,8.5))
show(pd)
}

\keyword{classes}


