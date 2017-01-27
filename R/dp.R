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
                nr=as.integer(nrow(M)), ncol=as.integer(ncol(M)),
          match=as.integer(match), nmatch=as.integer(0),
          PACKAGE="flagme"
          ) 
  #list(D=matrix(out$D,ncol=ncol(D)),phi=matrix(out$phi,ncol=ncol(D)),match=1+matrix(out$match,ncol=2)[out$nmatch:1,])
  list(match=1+matrix(out$match,ncol=2)[out$nmatch:1,])
}


## this two is to allow to use both the ndp
## functions and the matching functions from MR and RR

## .dynRT <- function(S){
##     options(warn=-1)
##     ## S similarity matrix
##     ## move trought S to find the highest score
##     S <- round(S, digits=3)
##     trace <- c()
##     for(i in 1:nrow(S)){
##         trace[i] <- which(S[i,] == max(S[i,]))    
##     }
##     trace.mtx <- matrix(NA, nrow=nrow(S), ncol=2)  
##     trace.mtx[,1] <- 1:nrow(S)
##     trace[which(trace == 1)] <- NA
##     trace.mtx[,2] <- trace

##     ## cercare i composti matchati più di una volta e renderli univoci
##     idx <- which(table(trace.mtx[,2]) > 1) # colonne della matrice S
##     names.idx <- as.numeric(names(idx)) # i doppioni
##     position <- c() # i doppi con il match più alto
##     for(k in 1:length(names.idx)){
##         position[k] <- which(S[,names.idx[k]] == max(S[,names.idx[k]]))
##     }
##     tutti <- which(trace %in% names.idx)
##     trace.mtx[,2][tutti[! tutti %in% position]] <- NA
        
##     return(list(match=trace.mtx))
    
##     options(warn=0)
## }

## .dynRTmin <- function(S){
##     options(warn=-1)
##     ## S similarity matrix
##     ## move trought S to find the highest score
##     S <- round(S, digits=3)
##     trace <- c()
##     for(i in 1:nrow(S)){
##         trace[i] <- which(S[i,] == min(S[i,]))##
##     }
##     trace.mtx <- matrix(NA, nrow=nrow(S), ncol=2)  
##     trace.mtx[,1] <- 1:nrow(S)
##     trace[which(trace == 1)] <- NA
##     trace.mtx[,2] <- trace

##     ## cercare i composti matchati più di una volta e renderli univoci
##     idx <- which(table(trace.mtx[,2]) > 1) # colonne della matrice S
##     names.idx <- as.numeric(names(idx)) # i doppioni
##     position <- c() # i doppi con il match più alto
##     for(k in 1:length(names.idx)){
##         position[k] <- which(S[,names.idx[k]] == min(S[,names.idx[k]]))
##     }
##     tutti <- which(trace %in% names.idx)
##     trace.mtx[,2][tutti[! tutti %in% position]] <- NA
        
##     return(list(match=trace.mtx))
    
##     options(warn=0)
## }



dynRT <- function(S){
    options(warn=-1)
    ## S similarity matrix
    ## move trought S to find the highest score
    S <- round(S, digits=3)
    filling <- which(table(S) > 200)
    filling.names <- as.numeric(names(filling))
    trace <- c()
    ## this boolean workaround is to allow to use both the ndp
    ## function and the matching function from MR and RR 
    if(filling.names == 1)
    {
        for(i in 1:nrow(S)){
            trace[i] <- which(S[i,] == min(S[i,])) # MR
        }
    }
    if(filling.names == 0)
    {
        for(i in 1:nrow(S)){
            trace[i] <- which(S[i,] == max(S[i,])) # RR
        }
    }
    
    trace.mtx <- matrix(NA, nrow=nrow(S), ncol=2)  
    trace.mtx[,1] <- 1:nrow(S)
    trace[which(trace == 1)] <- NA
    trace.mtx[,2] <- trace

    ## cercare i composti matchati più di una volta e renderli univoci
    idx <- which(table(trace.mtx[,2]) > 1) # colonne della matrice S
    names.idx <- as.numeric(names(idx)) # i doppioni
    position <- c() # i doppi con il match più alto

    if(filling.names == 1)
    {
        for(k in 1:length(names.idx)){
            position[k] <- which(S[,names.idx[k]] == min(S[,names.idx[k]]))
        } 
    }
    if(filling.names == 0)
    {
        for(k in 1:length(names.idx)){
            position[k] <- which(S[,names.idx[k]] == max(S[,names.idx[k]]))
        }
    }
 
    tutti <- which(trace %in% names.idx)
    trace.mtx[,2][tutti[! tutti %in% position]] <- NA
        
    return(list(match=trace.mtx))
    
    options(warn=0)
}
