\name{addXCMSPeaks}
\alias{addXCMSPeaks}

\title{Add xcms/CAMERA peak detection results}

\description{Reads the raw data using xcms, group each extracted ion
  according to their retention time using CAMERA and attaches them to an
  already created \code{peaksDataset} object} 

\usage{addXCMSPeaks(files, object, peakPicking=c('cwt','mF'), ...)} 

\arguments{
  \item{files}{character vector of same length as \code{object@rawdata}
    (user ensures the order matches)} 
  \item{object}{a \code{peaksDataset} object.}
  \item{peakPicking}{Methods to use for peak detection. See details.}
  \item{...}{arguments passed on to \code{xcmsSet} and \code{annotate}} 
}

\details{
  Repeated calls to xcmsSet and annotate to perform peak-picking and
  deconvolution. The peak detection results are added to the original
  \code{peaksDataset} object. Two peak detection alorithms are
  available: continuous wavelet transform (peakPicking=c('cwt')) and the
  matched filter approach (peakPicking=c('mF')) described by Smith et al
  (2006). For further information consult the xcms package manual.
}

\value{
  \code{peaksDataset} object
}

\author{Riccardo Romoli \email{riccardo.romoli@unifi.it}}

%% \references{
%%   Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data 
%%   \emph{PhD dissertation} University of Melbourne.
%% }

\seealso{
  \code{\link{peaksDataset}}
  \code{\link{findPeaks.matchedFilter}}
  \code{\link{findPeaks.centWave}}
  \code{\link{xcmsRaw-class}}
}

\examples{
# need access to CDF (raw data)
require(gcspikelite)
gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")

# full paths to file names
cdfFiles <- dir(gcmsPath, "CDF", full=TRUE)

# create a 'peaksDataset' object and add XCMS peaks to it
pd <- peaksDataset(cdfFiles[1], mz=seq(50,550), rtrange=c(7.5,8.5))
pd <- addXCMSPeaks(cdfFiles[1], pd, peakPicking=c('mF'),
                   snthresh=3, fwhm=4, step=1, steps=2, mzdiff=0.5)
}

\keyword{manip}


