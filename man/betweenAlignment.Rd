\name{betweenAlignment}

\alias{betweenAlignment}
\alias{betweenAlignment-class}
\alias{betweenAlignment-show}
\alias{show,betweenAlignment-method}

\title{Data Structure for "between" alignment of many GCMS samples}

\description{This function creates a "between" alignment (i.e. comparing merged peaks)}

\usage{betweenAlignment(pD,cAList,pAList,impList,filterMin=3,gap=0.7,D=10,usePeaks=TRUE,df=30,verbose=TRUE)}

\arguments{

\item{pD}{a \code{peaksDataset} object}

\item{cAList}{\code{list} of \code{clusterAlignment} objects, one for each experimental group}

\item{pAList}{\code{list} of \code{progressiveAlignment} objects, one for each experimental group}

\item{impList}{\code{list} of imputation lists}

\item{filterMin}{minimum number of peaks within a merged peak to be kept in the analysis}

\item{gap}{gap parameter}

\item{D}{retention time penalty parameter}

\item{usePeaks}{logical, whether to use peaks (if \code{TRUE}) or the full 2D profile alignment (if \code{FALSE})}

\item{df}{distance from diagonal to calculate similarity}

\item{verbose}{logical, whether to print information}

}


\details{

\code{betweenAlignment} objects gives the data structure which stores the result of an alignment across several "pseudo" datasets.  These pseudo datasets are constructed by merging the "within" alignments.

}

\value{

\code{betweenAlignment} object

}


\author{Mark Robinson}


\references{

Mark D Robinson (2008).  Methods for the analysis of gas chromatography - mass spectrometry data \emph{PhD dissertation} University of Melbourne.

}


\seealso{
\code{\link{multipleAlignment}}
}

\examples{
require(gcspikelite)
# see 'multipleAlignment'
}

\keyword{classes}


