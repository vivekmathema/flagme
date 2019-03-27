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


## RR ##
## retention time penalized normDotProd
ndpRT <- function(s1, s2, t1, t2, D){
    
    Normalize <- function(j){
        n <- apply(j, 2, function(k){
            m <- k[which.max(k)]
            norm <- k/m*100
        })
        return(n)
    }
    
    scoring <- function(s1, s2, t1, t2, D){
        angle <- function(s1, s2){
            theta <- acos(sum(s1*s2) / (sqrt(sum(s1 * s1)) * sqrt(sum(s2 * s2))))
            theta <- 1-theta
            if(theta < 0)
            {
                theta <- 0
            }       
            return(theta) 
        }
        
        rtPen <- function(t1, t2, D){
            ## D espresso in secondi
            t1 <- t1/60 # trasformo in secondi
            t2 <- t2/60 # trasformo in secondi
            srt <- exp(-(((t1-t2)^2) / D^2)) # da articolo MR, modificato
            # era 2*D^2
            return(srt)
        }
        
        score <- angle(s1, s2) * rtPen(t1, t2, D)
        return(score)
    }
    
    s1 <- Normalize(s1)
    s2 <- Normalize(s2)
    
    res <- matrix(0, nrow=ncol(s1), ncol=ncol(s2))
    for(i in 1:ncol(s1)){
        for(j in 1:ncol(s2)){
            res[i,j] <- scoring(s1[,i], s2[,j], t1[i], t2[j], D=D)
        }
    }    
    return(res)
}


## correlation Alignment
corPrt <- function(d1, d2, t1, t2, D, penality=0.2){
    D <- as.numeric(D) # time window in second
    pn <- as.numeric(penality)# penality if out of time window
    pearson <- function(x,y){
        size <- length(x)
        cfun <- .C("pearson", size=as.integer(size), x=as.double(x),
                   y=as.double(y), result=double(1), PACKAGE='flagme')
        return(cfun[["result"]])
    }
    Normalize <- function(j){
        n <- apply(j, 2, function(k){
            m <- k[which.max(k)]
            norm <- k/m*100
        })
    }
    Rank <- function(u) {
        if (length(u) == 0L) 
            u
        else if (is.matrix(u)) {
            if (nrow(u) > 1L) 
                apply(u, 2L, rank, na.last="keep")
            else row(u)
        }
        else rank(u, na.last="keep")
    }
    #
        x <- Normalize(d1)
        y <- Normalize(d2)
    
    ## method <- c("pearson", "kendall", "spearman")
    ncx <- ncol(x)
    ncy <- ncol(y)
    r <- matrix(0, nrow=ncx, ncol=ncy)
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncy)) {
            x2 <- x[, i]
            y2 <- y[, j]
            ok <- complete.cases(x2, y2)
            x2 <- rank(x2[ok])
            y2 <- rank(y2[ok])
            ## insert rt penality in seconds
            rtDiff <- t1[i]*60 - t2[j]*60 # retention time in seconds
            rtDiff <- abs(rtDiff)
            r[i, j] <- if (any(ok))
                           if(rtDiff <= D)
                               pearson(x2, y2)
                           else 
                               pearson(x2, y2) - pn
                       else 0
        }
    }
    
    r <- apply(r, MARGIN=c(1,2), function(x){
        if(x < 0.2)
        {
            x <- 0
        }
        else
        {
            x <- x
        }
    })
    rownames(r) <- colnames(x)
    colnames(r) <- colnames(y)
    return(r)
}
