#' Fits a robust linear model (RLM) for one metabolite
#' 
#' Using \code{rlm} from MASS, this procedure fits a linear model using all the
#' fragments
#' 
#' 
#' Fits a robust linear model.
#' 
#' @param u a metabolite unit (list object with vectors \code{mz} and \code{rt}
#' for m/z and retention times, respectively and a \code{data} element giving
#' the fragmentxsample intensitity matrix)
#' @param maxit maximum number of iterations (default: 5)
#' @param mzEffect logical, whether to fit m/z effect (default: \code{TRUE})
#' @param cls class variable
#' @param fitSample whether to fit individual samples (alternative is fit by
#' group)
#' @param fitOrCoef whether to return a vector of coefficients (default:
#' "coef"), or an \code{rlm} object ("fit")
#' @param TRANSFORM function to transform the raw data to before fitting
#' (default: \code{log2})
#' @return
#' 
#' \code{list} giving elements of \code{fragment} and \code{sample}
#' coefficients (if \code{fitOrCoef="coef"}) or a \code{list} of elements from
#' the fitting process (if \code{fitOrCoef="fit"})
#' @author Mark Robinson
#' @seealso \code{\link{peaksAlignment}}, \code{\link{clusterAlignment}}
#' @references
#' 
#' Mark D Robinson (2008).  Methods for the analysis of gas chromatography -
#' mass spectrometry data \emph{PhD dissertation} University of Melbourne.
#' @keywords manip
#' @examples
#' 
#' require(gcspikelite)
#' 
#' # paths and files
#' gcmsPath<-paste(find.package("gcspikelite"),"data",sep="/")
#' cdfFiles<-dir(gcmsPath,"CDF",full=TRUE)
#' eluFiles<-dir(gcmsPath,"ELU",full=TRUE)
#' 
#' # read data, peak detection results
#' pd<-peaksDataset(cdfFiles[1:2],mz=seq(50,550),rtrange=c(7.5,8.5))
#' pd<-addAMDISPeaks(pd,eluFiles[1:2])
#' 
#' # pairwise alignment using all scans
#' fullca<-clusterAlignment(pd, usePeaks = FALSE, df = 100)
#' 
#' # calculate retention time shifts
#' timedf<-calcTimeDiffs(pd, fullca)
#'
#' @importFrom stats contr.sum model.matrix
#' @importFrom MASS rlm
#' @export rmaFitUnit
rmaFitUnit <- function(u,maxit=5,mzEffect=TRUE,cls=NULL,fitSample=TRUE,fitOrCoef=c("coef","fit"),TRANSFORM=log2) {
  d<-u$data
  fitOrCoef=match.arg(fitOrCoef)
  k<-rowSums(d==0 | is.na(d))==0  # cannot have rows with 0s and then take logs
  #if( any(is.na(d)) ) {
  #  cat('Unit has NAs.  Please check\n')
  #  return(NULL)
  #}
  y <- as.vector(TRANSFORM(d[k,]))
  probe <- factor(rep(u$mz[k],ncol(d)))
  sample <- factor(rep(names(u$rt),each=nrow(d[k,])))
  clsv <- factor(rep(cls,each=nrow(d[k,])))
  if( !mzEffect ) {
    if (fitSample)
      v<-model.matrix(~-1+sample+probe,contrasts=list(probe=contr.sum))
    else
      v<-model.matrix(~-1+clsv+probe,contrasts=list(probe=contr.sum))
  } else {
    if ( is.null(cls) )
      stop("'cls' needs to be specified in order to run this model.")
    mz<-rep(u$mz[k],ncol(d))
    if (fitSample) {
          vv<-model.matrix(~-1+sample+probe,contrasts=list(probe=contr.sum))
          vv1<-model.matrix(~-1+mz:clsv,contrast=list(clsv=contr.sum))
          dd<- vv1 %*% attr(vv1,"contrasts")$clsv
          v<-cbind(vv,dd)
        } else {
          vv<-model.matrix(~-1+clsv+probe,contrasts=list(probe=contr.sum))
          vv1<-model.matrix(~-1+mz:clsv,contrast=list(clsv=contr.sum))
          dd<- vv1 %*% attr(vv1,"contrasts")$clsv
          v<-cbind(vv,dd)
        }
    #else
    #  v<-model.matrix(~-1+clsv+probe+mz:clsv,contrasts=list(probe=contr.sum,clsv=contr.sum))
    #v<-v[,-ncol(v)]  # make design matrix full rank
  }
  fit<-rlm(y~-1+v,maxit=maxit)
  cat(".")
  if (fitOrCoef=="coef") {
    mzeff <- NULL
    sample<-fit$coef[1:ncol(d)]
        fragment<-fit$coef[(ncol(d)+1):(ncol(d)+nrow(d[k,])-1)]
        fragment<-c(fragment,-sum(fragment))
        if(mzEffect) {
		  mzeffects <- c(fit$coef[(ncol(d)+nrow(d[k,])):length(fit$coef)],0)
          mzeff<-mzeffects*mean(u$mz[k])
          sample<-sample+mzeff[sort(factor(cls))]
        }
        list(sample=sample,fragment=fragment,mzeff=mzeffects)
  } else {
    # save a little space for stuff that probably don't need
    #fit$qr<-NULL
        #fit$model<-NULL
        #fit$weights<-NULL
        #fit$effects<-NULL
        #fit$fitted.values<-NULL
    list(fit=fit,v=v,y=y)
  }
}
