\name{ndpRT}
\alias{ndpRT}

\title{Retention Time Penalized Correlation}

\description{This function calculates the similarity of all pairs of
  peaks from 2 samples, using the spectra similarity and the rretention
  time differencies}

\usage{corPrt(d1, d2, t1, t2, D=D, penality=0.2)}

\arguments{
\item{d1}{data matrix for sample 1}
\item{d2}{data matrix for sample 2}
\item{t1}{vector of retention times for sample 1}
\item{t2}{vector of retention times for sample 2}
\item{D}{retention time window for the matching}
\item{penality}{penalization applied to the matching between two mass spectra if \code{(t1-t2)>D}}
}

\details{
  Computes the Pearson carrelation between every pair of peak vectors
  in the retention time window (\code{D})and returns the similarity matrix.
}

\value{
  matrix of similarities
}

\author{Riccardo Romoli}

%% \references{
%%   Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data 
%%   \emph{PhD dissertation} University of Melbourne.
%% }

\seealso{
  \code{\link{peaksAlignment}}
}

\examples{
## Not Run
require(gcspikelite)
gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
cdfFiles <- dir(gcmsPath,"CDF", full=TRUE)
## read data, peak detection results
pd <- peaksDataset(cdfFiles[1:3], mz=seq(50,550), rtrange=c(7.5,10.5))
pd <- addXCMSPeaks(files=cdfFiles[1:3], object=pd, peakPicking=c('mF'),
                   snthresh=3, fwhm=10,  step=0.1, steps=2, mzdiff=0.5,
                   sleep=0)
## review peak picking
plot(pd, rtrange=c(7.5, 10.5), runs=c(1:3))

r <- corPrt(pd@peaksdata[[1]], pd@peaksdata[[2]],
           pd@peaksrt[[1]], pd@peaksrt[[2]], D=50, penality=0.2)
## End (Not Run)
}

\keyword{manip}
