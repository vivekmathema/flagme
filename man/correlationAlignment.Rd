\name{correlatioAlignment}

%% \alias{correlationAlignment-class}
%% \alias{correlationAlignment-show}
\alias{correlationAlignment}
%% \alias{show, correlationAlignment-method}
%% \alias{show, correlationAlignment, ANY-method}

\title{Function for correlation-based alignment strategy of multiple GCMS samples}

\description{Performs a correlation-based peak alignment of multiple
  GCMS peak lists; this function is able to align multiple samples, by a
  center-star strategy.}

\usage{
correlationAlignment(object, thr=0.85, D=20, penality=0.2, normalize=TRUE, minFilter=1)
} 

\arguments{
  \item{object}{a \code{peaksDataset} object}
  \item{thr}{correlation trashold from 0 (min) to 1 (max)}
  \item{D}{retention time window in seconds}
  \item{penality}{the penality inflicted to a match between two
    peaks when the retention time difference exceed the parameter D}
  \item{normalize}{logical, whether to use normalized-to-100 peaks
  intensity or as such}
  \item{minFilter}{if a feature is matched in a number of samples less
    than \code{minFilter}, this feature is trashed. The value of minFilter
    must be smaller than the number of samples}
}
  
\details{
  The correlation-based peak alignment for multiple GCMS
  peak lists uses a center-star tecnique to the alignment of the
  peaks. The combination of the \code{D} and \code{penality} parameters
  allow the users to force the algorithm to match the peaks close to the
  reference. The \code{thr} parameter control the matching factor.
}

\value{
  \code{correlationAlignment} object
}

\author{Riccardo Romoli \email{riccardo.romoli@unifi.it}}

\seealso{
  \code{\link{peaksDataset}}, \code{\link{addXCMSPeaks}},
    \code{\link{correlationAlignment-class}
  }
}

\examples{
require(gcspikelite)

# paths and files
gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
cdfFiles <- dir(gcmsPath,"CDF",full=TRUE)

# read data, peak detection results
pd <- peaksDataset(cdfFiles[1:2], mz=seq(50,550), rtrange=c(7.5,8.5))
pd <- addXCMSPeaks(files=cdfFiles[1:2], object=pd, peakPicking=c('mF'),
                   snthresh=3, fwhm=4, step=1, steps=2, mzdiff=0.5)
mp <- correlationAlignment(object=pd, thr=0.85, D=20, penality=0.2,
                           normalize=TRUE, minFilter=1)
}

\keyword{classes}