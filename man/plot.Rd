\name{plot.peaksDataset}

\alias{plot}
\alias{.plotpA}
\alias{.plotpD}
\alias{.plotcA}

\title{Plotting functions for GCMS data objects}

\description{Store the raw data and optionally, information regarding
  signal peaks for a number of GCMS runs} 

\usage{
.plotpD(object, runs=1:length(object@rawdata),
       mzind=1:nrow(object@rawdata[[1]]), mind=NULL,
       plotSampleLabels=TRUE, calcGlobalMax=FALSE, peakCex = 0.8,
       plotPeaks=TRUE, plotPeakBoundaries=FALSE, plotPeakLabels=FALSE,
       plotMergedPeakLabels=TRUE, mlwd=3,usePeaks=TRUE,
       plotAcrossRuns=FALSE, overlap=F, rtrange=NULL, cols=NULL, thin=1,
       max.near=median(object@rawrt[[1]]), how.near=50, scale.up=1, ...)
			   
.plotpA(object, xlab="Peaks - run 1", ylab="Peaks - run 2",
       plotMatches=TRUE, matchPch=19, matchLwd=3,
       matchCex=.5, matchCol="black", col=colorpanel(50, "white", "green", "navyblue"),
%% col=colorpanel(50,"black","blue","white"),
       breaks=seq(0, 1, length=51), ...)
			   
.plotcA(object, alignment=1, ...)
}

\arguments{
  \item{object}{a \code{peaksDataset}, \code{peaksAlignment} or \code{clusterAlignment} object.}
  \item{runs}{for \code{peaksDataset} only: set of run indices to plot}
  \item{mzind}{for \code{peaksDataset} only: set of mass-to-charge indices to sum over (default, all)}
  \item{mind}{for \code{peaksDataset} only: matrix of aligned indices}
  \item{plotSampleLabels}{for \code{peaksDataset} only: logical, whether to display sample labels}
  \item{calcGlobalMax}{for \code{peaksDataset} only: logical, whether to calculate an overall maximum for scaling}
  \item{peakCex}{character expansion factor for peak labels}
  \item{plotPeaks}{for \code{peaksDataset} only: logical, whether to plot hashes for each peak}
  \item{plotPeakBoundaries}{for \code{peaksDataset} only: logical, whether to display peak boundaries}
  \item{plotPeakLabels}{for \code{peaksDataset} only: logical, whether to display peak labels}
  \item{plotMergedPeakLabels}{for \code{peaksDataset} only: logical, whether to display 'merged' peak labels}
  \item{mlwd}{for \code{peaksDataset} only: line width of lines indicating the alignment}
  \item{usePeaks}{for \code{peaksDataset} only: logical, whether to plot alignment of peaks (otherwise, scans)}
  \item{plotAcrossRuns}{for \code{peaksDataset} only: logical, whether to plot across peaks when unmatched peak is given}
  \item{overlap}{for \code{peaksDataset} only: logical, whether to plot TIC/XICs overlapping}
  \item{rtrange}{for \code{peaksDataset} only: vector of length 2 giving start and end of the X-axis}
  \item{cols}{for \code{peaksDataset} only: vector of colours (same length as the length of runs)}
  \item{thin}{for \code{peaksDataset} only: when \code{usePeaks=FALSE}, plot the alignment lines every \code{thin} values}
  \item{max.near}{for \code{peaksDataset} only: where to look for maximum}
  \item{how.near}{for \code{peaksDataset} only: how far away from \code{max.near} to look}
  \item{scale.up}{for \code{peaksDataset} only: a constant factor to scale the TICs}
  \item{plotMatches}{for \code{peaksDataset} only: logical, whether to plot matches}
  \item{xlab}{for \code{peaksAlignment} and \code{clusterAlignment} only: x-axis label}
  \item{ylab}{for \code{peaksAlignment} and \code{clusterAlignment} only: y-axis label}
  \item{matchPch}{for \code{peaksAlignment} and \code{clusterAlignment} only: match plotting character}
  \item{matchLwd}{for \code{peaksAlignment} and \code{clusterAlignment} only: match line width}
  \item{matchCex}{for \code{peaksAlignment} and \code{clusterAlignment} only: match character expansion factor}
  \item{matchCol}{for \code{peaksAlignment} and \code{clusterAlignment} only: match colour}
  \item{col}{for \code{peaksAlignment} and \code{clusterAlignment} only: vector of colours for colourscale}
  \item{breaks}{for \code{peaksAlignment} and \code{clusterAlignment} only: vector of breaks for colourscale}
  \item{alignment}{for \code{peaksAlignment} and \code{clusterAlignment} only: the set of alignments to plot}
  \item{...}{further arguments passed to the \code{plot} or \code{image} command}
}

\details{
  For \code{peakDataset} objects, each TIC is scale to the maximum value
  (as specified by the \code{how.near} and \code{max.near} values).  The
  many parameters gives considerable flexibility of how the TICs can be
  visualized. 
  
  For \code{peakAlignment} objects, the similarity matrix is plotted and
  optionally, the set of matching peaks.  \code{clusterAlignment}
  objects are just a collection of all pairwise \code{peakAlignment}
  objects. 
}

\author{Mark Robinson}

\references{
Mark D Robinson (2008).  Methods for the analysis of gas chromatography
- mass spectrometry data \emph{PhD dissertation} University of
Melbourne. 
}

\seealso{
  \code{\link{plotImage}}, \code{\link{peaksDataset}}
}

\examples{
require(gcspikelite)

## paths and files
gcmsPath <- paste(find.package("gcspikelite"), "data", sep="/")
cdfFiles <- dir(gcmsPath, "CDF", full=TRUE)
eluFiles <- dir(gcmsPath, "ELU", full=TRUE)

## read data
pd <- peaksDataset(cdfFiles[1:3], mz=seq(50,550), rtrange=c(7.5,8.5))

## image plot
plot(pd, rtrange=c(7.5,8.5), plotPeaks=TRUE, plotPeakLabels=TRUE)
}

\keyword{classes}


