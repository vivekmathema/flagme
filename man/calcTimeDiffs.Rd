\name{calcTimeDiffs}
\alias{calcTimeDiffs}

\title{Calculate retention time shifts from profile alignments}

\description{This function takes the set of all pairwise profile alignments and use these to estimate retention time shifts between each pair of samples.  These will then be used to normalize the retention time penalty of the signal peak alignment.}

\usage{calcTimeDiffs(pd,ca.full,verbose=TRUE)}

\arguments{

\item{pd}{a \code{peaksDataset} object}

\item{ca.full}{a \code{clusterAlignment} object, fit with }

\item{verbose}{logical, whether to print out information}

}


\details{

Using the set of profile alignments, 

}

\value{

\code{list} of same length as \code{ca.full@alignments} with the matrices giving the retention time penalties.

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data 
\emph{PhD dissertation} University of Melbourne.

}


\seealso{
\code{\link{peaksAlignment}}, \code{\link{clusterAlignment}}
}

\examples{
require(gcspikelite)

# paths and files
gcmsPath<-paste(.find.package("gcspikelite"),"data",sep="/")
cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
eluFiles<-dir(gcmsPath,"ELU",full=TRUE)

# read data, peak detection results
pd<-peaksDataset(cdfFiles[1:2],mz=seq(50,550),rtrange=c(7.5,8.5))
pd<-addAMDISPeaks(pd,eluFiles[1:2])

# pairwise alignment using all scans
fullca<-clusterAlignment(pd, usePeaks = FALSE, df = 100)

# calculate retention time shifts
timedf<-calcTimeDiffs(pd, fullca)
}

\keyword{manip}


