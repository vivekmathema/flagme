\name{dp}
\alias{dp}

\title{Dynamic programming algorithm, given a similarity matrix}

\description{This function calls C code for a bare-bones dynamic programming algorithm, finding the best cost path through a similarity matrix.}

\usage{dp(M,gap=.5,big=10000000000,verbose=FALSE)}

\arguments{

\item{M}{similarity matrix}

\item{gap}{penalty for gaps}

\item{big}{large value used for matrix margins}

\item{verbose}{logical, whether to print out information}

}


\details{

This is a pretty standard implementation of a bare-bones dynamic programming algorithm, with a single gap parameter and allowing only simple jumps through the matrix (up, right or diagonal).

}

\value{

\code{list} with element \code{match} with the set of pairwise matches.

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data 
\emph{PhD dissertation} University of Melbourne.

}


\seealso{
\code{\link{normDotProduct}}
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

# similarity matrix
r<-normDotProduct(pd@peaksdata[[1]],pd@peaksdata[[2]])

# dynamic-programming-based matching of peaks
v<-dp(r,gap=.5)
}

\keyword{manip}


