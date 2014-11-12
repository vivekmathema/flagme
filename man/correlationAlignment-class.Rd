\name{correlationAlignment-class}
\docType{class}

\alias{correlationAlignment-class}
\alias{correlationAlignment-show}
\alias{show,correlationAlignment-method}
\alias{show,correlationAlignment-method, ANY-methods}
\alias{align,correlatioAligment}

\title{\code{"correlationAlignment-class"}}
\description{A class containing the index of the aligned chromatographic peaks.}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("correlationAlignment", ...)}. 
  The object created contains both the results of the alignment procedure
  and the file used and as center-star.
}

\section{Slots}{
  \describe{
    \item{\code{Alignment}:}{Object of class \code{"align"} contain the
      matrix of the aligned features. The rows represent the different
      peaks while the columns represent the files. The values of the
      matrix refers to \code{peaksind} slot of the \code{peaksDataSet}
      object. }
    \item{\code{Center}:}{Object of class \code{"character"} contain the
      file used as a center-star.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "correlationAlignment")}}
	 }
}

\author{Riccardo Romoli {\email{riccardo.romoli@unifi.it}}}

\seealso{%% \code{\linkS4class{addXCMSPeaks}} for detecting peaks and
  \code{\linkS4class{correlationAlignment}} for further information about the
  alignment function.
}

\examples{
showClass("correlationAlignment")
}

\keyword{classes}
