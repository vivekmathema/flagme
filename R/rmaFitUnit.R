.rmaFitPeak <- function(u,maxit=5,mzEffect=TRUE,cls=NULL,fitSample=TRUE,fitOrCoef=c("coef","fit"),TRANSFORM=log2) {
  require(MASS)
  d<-u$data
  fitOrCoef=match.arg(fitOrCoef)
  k<-rowSums(d==0 | is.na(d))==0  # cannot have rows with 0s and then take logs
  #if( any(is.na(d)) ) {
  #  cat('Unit has NAs.  Please check\n')
  #  return(NULL)
  #}
  y<-as.vector(TRANSFORM(d[k,]))
  probe<-factor(rep(u$mz[k],ncol(d)))
  sample<-factor(rep(names(u$rt),each=nrow(d[k,])))
  clsv<-factor(rep(cls,each=nrow(d[k,])))
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
