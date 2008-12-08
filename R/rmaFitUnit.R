.rmaFitPeak <- function(u,PsiCode=0,PsiK=1.345) {
  d<-u$data
  d<-d[rowSums(is.na(d) | d==0)==0,]
  y<-log(d,base=2)
  beta<-rep(NA,ncol(y))
  keep<-colSums(is.na(y))==0
  y<-y[,keep]
  fit<-.Call("R_rlm_rma_default_model",y,PsiCode,PsiK,PACKAGE="preprocessCore")
  J <- ncol(y)
  est <- fit$Estimates
  #se <- fit$StdErrors
  beta.tmp <- est[1:J]
  alpha <- est[(J + 1):length(est)]
  alpha[length(alpha)] <- -sum(alpha[1:(length(alpha)-1)])
  beta[keep]<-beta.tmp
  list(alpha=alpha,beta=beta,res=fit$Residuals,keep=keep)
}

