normDotProduct <- function (x1, x2, t1=NULL, t2=NULL, df=max(ncol(x1), ncol(x2)), 
                               D=1e+05, timedf=NULL, verbose=FALSE){
    if(is.null(t1)) 
        t1 <- 1:ncol(x1)
    if(is.null(t2)) 
        t2 <- 1:ncol(x2)
    if(nrow(x1) != nrow(x2) | length(t1) != ncol(x1) | length(t2) != 
       ncol(x2)){
        stop("One of these is not true: nrow(x1)=nrow(x2), length(t1)=ncol(x1), length(t2)=ncol(x2).")
    }
    score <- matrix(0, nrow=ncol(x1), ncol=ncol(x2))
    if(length(D) == 1 & is.null(timedf)){
        out <- .C("cos_ndp_lowmem", score=as.double(score), 
                  nr=as.integer(nrow(x1)), nc1=as.integer(ncol(x1)), 
                  nc2=as.integer(ncol(x2)), x1=as.double(x1), x2=as.double(x2), 
                  t1=as.double(t1), t2 = as.double(t2), D = as.double(D), 
                  df=as.integer(df), PACKAGE="flagme")
    }else{
        if(length(D) == 1) 
            D <- matrix(D, nrow=ncol(x1), ncol=ncol(x2))
        if(ncol(x1) != nrow(D) | ncol(x2) != ncol(D)) 
            stop("D must have dimensions nrow=ncol(x1) ncol=ncol(x2) or be scalar.")
        if(is.null(timedf)) 
            timedf <- matrix(0, nrow=ncol(x1), ncol=ncol(x2))
        if(ncol(x1) != nrow(timedf) | ncol(x2) != ncol(timedf)) 
            stop("'timedf' must have dimensions nrow=ncol(x1) ncol=ncol(x2) or be set to NULL.")
        out <- .C("cos_ndp_himem", score=as.double(score), 
                  nr=as.integer(nrow(x1)), nc1=as.integer(ncol(x1)), 
                  nc2=as.integer(ncol(x2)), x1=as.double(x1), x2=as.double(x2), 
                  D=as.double(D), df=as.integer(df), timedf=as.double(timedf), 
                  PACKAGE="flagme")
    }
    NDP <- matrix(1 - out$score, ncol=ncol(x2))
    NDP[is.nan(NDP)] <- 0 ## remove NaN    
    return(NDP)
}

## normDotProduct<-function(x1,x2,t1=NULL,t2=NULL,df=max(ncol(x1),ncol(x2)),D=100000,timedf=NULL,verbose=FALSE){
##   if(is.null(t1)) t1<-1:ncol(x1)
##   if(is.null(t2)) t2<-1:ncol(x2)
  
##   if( nrow(x1) != nrow(x2) | length(t1)!=ncol(x1) | length(t2)!=ncol(x2) ) {
##     stop("One of these is not true: nrow(x1)=nrow(x2), length(t1)=ncol(x1), length(t2)=ncol(x2).")
##   }

##   score<-matrix(0,nrow=ncol(x1),ncol=ncol(x2))
  
##   if( length(D)==1 & is.null(timedf) ) {
##     out<-.C("cos_ndp_lowmem", score=as.double(score), nr=as.integer(nrow(x1)), nc1=as.integer(ncol(x1)),
##                      nc2=as.integer(ncol(x2)), x1=as.double(x1), x2=as.double(x2), t1=as.double(t1),
##                      t2=as.double(t2),D=as.double(D),df=as.integer(df), PACKAGE="flagme")
##   } else {
##     if (length(D)==1)
## 	  D<-matrix(D,nrow=ncol(x1),ncol=ncol(x2))
##     if( ncol(x1) != nrow(D) | ncol(x2)!=ncol(D) )
##       stop("D must have dimensions nrow=ncol(x1) ncol=ncol(x2) or be scalar.")
## 	if( is.null(timedf) )
## 	  timedf<-matrix(0,nrow=ncol(x1),ncol=ncol(x2))
##     if( ncol(x1) != nrow(timedf) | ncol(x2)!=ncol(timedf) )
##       stop("'timedf' must have dimensions nrow=ncol(x1) ncol=ncol(x2) or be set to NULL.")
##     out<-.C("cos_ndp_himem", score=as.double(score), nr=as.integer(nrow(x1)), nc1=as.integer(ncol(x1)),
##                      nc2=as.integer(ncol(x2)), x1=as.double(x1), x2=as.double(x2),
## 					 D=as.double(D),df=as.integer(df),timedf=as.double(timedf), PACKAGE="flagme")    
##   }
##   matrix(1-out$score,ncol=ncol(x2))
## }
