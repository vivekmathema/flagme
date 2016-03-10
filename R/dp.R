dp<-function(M,gap=.5,big=10000000000,verbose=FALSE) {

  # setup score matrix
  bigr<-c(0,big*(1:ncol(M)))
  bigc<-c(big*(1:nrow(M)))
  D<-rbind(bigr,cbind(bigc,M)); #D[1,1]<-0
  #bigr<-c(0,gap*(1:ncol(M)))
  #bigc<-c(gap*(1:nrow(M)))
  #D<-rbind(bigr,cbind(bigc,M)); #D[1,1]<-0
  
  # setup traceback
  phi<-matrix(0,nrow=nrow(D),ncol=ncol(D))
  phi[1,]<-2; phi[,1]<-1; phi[1,1]<-3
  
  match<-matrix(-1,nrow=nrow(M)+ncol(M),ncol=2) # matrix to store matches

  out<-.C("dp", D=as.double(D), M=as.double(M), gap=as.double(gap), phi=as.integer(phi), 
                nr=as.integer(nrow(M)), ncol=as.integer(ncol(M)), match=as.integer(match), nmatch=as.integer(0), PACKAGE="flagme")
  #list(D=matrix(out$D,ncol=ncol(D)),phi=matrix(out$phi,ncol=ncol(D)),match=1+matrix(out$match,ncol=2)[out$nmatch:1,])
  list(match=1+matrix(out$match,ncol=2)[out$nmatch:1,])
}
