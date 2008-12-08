
.plotFragments<-function(peaks,pk=1,cols=NULL,ltys=NULL,TRANSFUN=log2,...) {
  if(is.null(cols))
    cols<-rep(c("red","grey40","green","blue"),each=ncol(peaks[[pk]]$data)/4)
  if(is.null(ltys))
     ltys<-rep(1,length(peaks))
  pr<-paste(round(range(peaks[[pk]]$rt,na.rm=TRUE),2),collapse="-")
  matplot(peaks[[pk]]$mz,TRANSFUN(peaks[[pk]]$data),type="l",col=cols,lty=ltys,xlab="m/z",ylab="Intensity",main=paste("index:",pk,"rt:",pr, "fragments:",nrow(peaks[[pk]]$data)),...)
}

