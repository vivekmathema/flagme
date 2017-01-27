\name{ndpRT}
\alias{ndpRT}

\title{Retention Time Penalized Normalized Dot Product}

\description{This function calculates the similarity of all pairs of
  peaks from 2 samples, using the spectra similarity and the rretention
  time differencies}

\usage{ndpRT(s1, s2, t1, t2, D)}

\arguments{
\item{s1}{data matrix for sample 1}
\item{s2}{data matrix for sample 2}
\item{t1}{vector of retention times for sample 1}
\item{t2}{vector of retention times for sample 2}
\item{D}{retention time window for the matching}
}

\details{
  Computes the normalized dot product between every pair of peak vectors
  in the retention time window (\code{D})and returns a similarity matrix.
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

                                        # read data, peak detection results
pd <- peaksDataset(cdfFiles[1:3], mz=seq(50,550), rtrange=c(7.5,10.5))
pd <- addXCMSPeaks(files=cdfFiles[1:3], object=pd, peakPicking=c('mF'),
                   snthresh=3, fwhm=10,  step=0.1, steps=2, mzdiff=0.5,
                   sleep=0)
## review peak picking
plot(pd, rtrange=c(7.5, 10.5), runs=c(1:3))

r <- ndpRT(pd@peaksdata[[1]], pd@peaksdata[[2]],
           pd@peaksrt[[1]], pd@peaksrt[[2]], D=50)
## End (Not Run)
}

\keyword{manip}
