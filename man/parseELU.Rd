\name{parseELU}
\alias{parseELU}

\title{Parser for ELU files}

\description{Reads ASCII ELU-format files from AMDIS (Automated Mass Spectral Deconvolution and Identification System)}

\usage{parseELU(f,min.pc=.01,mz=seq(50,550),rt.cut=.008,rtrange=NULL)}

\arguments{

\item{f}{ELU filename to read.}

\item{min.pc}{minimum percent of maximum intensity.}

\item{mz}{vector of mass-to-charge bins of raw data table.}

\item{rt.cut}{the difference in retention time, below which peaks are merged together.}

\item{rtrange}{retention time range to parse peaks from, can speed up parsing if only interested in a small region (must be \code{numeric} vector of length 2)}

}


\details{

\code{parseELU} will typically be called by \code{\link{addAMDISPeaks}}, not called directly.

Peaks that are detected within \code{rt.cut} are merged together.  This avoids peaks which are essentially
overlapping.

Fragments that are less than \code{min.pc} of the maximum intensity fragment are discarded.

}

\value{

\code{list} with components \code{peaks} (table of spectra -- rows are mass-to-charge and columns are the 
different detected peaks) and \code{tab} (table of features for each detection), according to what is stored in 
the ELU file.

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data 
\emph{PhD dissertation} University of Melbourne.

}


\seealso{
\code{\link{addAMDISPeaks}}
}

\examples{
require(gcspikelite)

# paths and files
gcmsPath<-paste(.find.package("gcspikelite"),"data",sep="/")
eluFiles<-dir(gcmsPath,"ELU",full=TRUE)

# parse ELU file
eluList<-parseELU(eluFiles[1])
}

\keyword{manip}


