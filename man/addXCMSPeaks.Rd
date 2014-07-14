\name{addXCMSPeaks}
\alias{addXCMSPeaks}

\title{Add xcms/CAMERA peak detection results}

\description{Reads the raw data using xcms, group each extracted ion
  according to their retention time using CAMERA and attaches them to an
  already created \code{peaksDataset} object} 

\usage{addXCMSPeaks(files, object, ...)} 

\arguments{
  \item{files}{character vector of same length as \code{object@rawdata}
    (user ensures the order matches)} 
  \item{object}{a \code{peaksDataset} object.} 
  \item{...}{arguments passed on to \code{xcmsSet} and \code{annotate}} 
}

\details{
  Repeated calls to xcmsSet and annotate to perform peak-picking and
  deconvolution. The peak detection results are added to the original
  \code{peaksDataset} object. Continuous wavelet transform (CWT) is used to
  locate chromatographic peaks on different scales.
}

\value{
  \code{peaksDataset} object
}

\author{Riccardo Romoli}

%% \references{
%%   Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data 
%%   \emph{PhD dissertation} University of Melbourne.
%% }

\seealso{
  \code{\link{peaksDataset}}
}

\examples{
# need access to CDF (raw data)
require(gcspikelite)
gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")

# full paths to file names
cdfFiles <- dir(gcmsPath, "CDF", full=TRUE)

# create a 'peaksDataset' object and add XCMS peaks to it
pd <- peaksDataset(cdfFiles[1], mz=seq(50,550), rtrange=c(7.5,8.5))
pd <- addXCMSPeaks(cdfFiles[1], pd)
}

\keyword{manip}