\name{peaksAlignment-class}

\alias{peaksAlignment-class}
\alias{peaksAlignment-show}
\alias{peaksAlignment-plot}
\alias{peaksAlignment}
\alias{show,peaksAlignment-method}
\alias{plot,peaksAlignment-method}
\alias{plot,peaksAlignment,ANY-method}

\title{Data Structure for pairwise alignment of 2 GCMS samples}

\description{Store the raw data and optionally, information regarding signal peaks for a number of GCMS runs}

\usage{peaksAlignment(d1,d2,t1,t2,gap=.5,D=1000,timedf=NULL,df=30,verbose=TRUE,usePeaks=TRUE,compress=TRUE)}

\arguments{

\item{d1}{matrix of MS intensities for 1st sample (if doing a peak alignment, this contains peak apexes/areas; if doing a profile alignment, this contains scan intensities.  Rows are m/z bins, columns are peaks/scans.}

\item{d2}{matrix of MS intensities for 2nd sample}

\item{t1}{vector of retention times for 1st sample}

\item{t2}{vector of retention times for 2nd sample}

\item{gap}{gap penalty for dynamic programming algorithm}

\item{D}{time penalty (on same scale as retention time differences, \code{t1} and \code{t2})}

\item{timedf}{list (length = the number of pairwise alignments) of matrices giving the expected time differences expected at each pair of peaks (used with \code{usePeaks}=\code{TRUE}.}

\item{df}{integer, how far from the diagonal to go to calculate the similarity of peaks.  Smaller value should run faster, but be careful not to choose too low.}

\item{verbose}{logical, whether to print out info.}

\item{usePeaks}{logical, \code{TRUE} uses \code{peakdata} list, \code{FALSE} uses \code{rawdata} list for computing similarity.}

\item{compress}{logical, whether to compress the similarity matrix into a sparse format.}

}


\details{

peaksAlignment is a hold-all data structure of the raw and peak detection data.

}

\value{

\code{peaksAlignment} object

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data \emph{PhD dissertation} University of Melbourne.

}


\seealso{
\code{\link{peaksDataset}}, \code{\link{clusterAlignment}}
}

\examples{

# see clusterAlignment, it calls peaksAlignment

}

\keyword{classes}


