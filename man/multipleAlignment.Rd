\name{multipleAlignment-class}

\alias{multipleAlignment-class}
\alias{multipleAlignment-show}
\alias{multipleAlignment}
\alias{show,multipleAlignment-method}

\title{Data Structure for multiple alignment of many GCMS samples}

\description{Store the raw data and optionally, information regarding signal peaks for a number of GCMS runs}

\usage{multipleAlignment(pd,group,bw.gap=0.8,wn.gap=0.6,bw.D=.20,wn.D=.05,filterMin=3,lite=FALSE,usePeaks=TRUE,df=50,verbose=TRUE,timeAdjust=FALSE,doImpute=FALSE)}

\arguments{

\item{pd}{a \code{peaksDataset} object}

\item{group}{factor variable of experiment groups, used to guide the alignment algorithm}

\item{bw.gap}{gap parameter for "between" alignments}

\item{wn.gap}{gap parameter for "within" alignments}

\item{bw.D}{distance penalty for "between" alignments}

\item{wn.D}{distance penalty for "within" alignments}

\item{filterMin}{minimum number of peaks within a merged peak to be kept in the analysis}

\item{lite}{logical, whether to keep "between" alignment details (default, \code{FALSE})}

\item{usePeaks}{logical, whether to use peaks (if \code{TRUE}) or the full 2D profile alignment (if \code{FALSE})}

\item{df}{distance from diagonal to calculate similarity}

\item{verbose}{logical, whether to print information}

\item{timeAdjust}{logical, whether to use the full 2D profile data to estimate retention time drifts (Note: time required)}

\item{doImpute}{logical, whether to impute the location of unmatched peaks}
}


\details{

multipleAlignment is the data structure giving the result of an alignment across several GCMS runs.

Multiple alignments are done progressively.  First, all samples with the same \code{tg$Group} label with be aligned (denoted a "within" alignment).  Second, each group will be summarized into a pseudo-data set, essentially a spectrum and retention time for each matched peak of the within-alignment.  Third, these "merged peaks" are aligned in the same progressive manner, here called a "between" alignment.

}

\value{

\code{multipleAlignment} object

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data \emph{PhD dissertation} University of Melbourne.

}


\seealso{
\code{\link{peaksDataset}}, \code{\link{betweenAlignment}}, \code{\link{progressiveAlignment}}
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

# multiple alignment
ma<-multipleAlignment(pd,c(1,1),wn.gap=0.5,wn.D=.05,bw.gap=0.6,bw.D=.2,usePeaks=TRUE,filterMin=1,df=50,verbose=TRUE)
}

\keyword{classes}


