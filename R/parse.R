#' Parser for ELU files
#' 
#' Reads ASCII ELU-format files from AMDIS (Automated Mass Spectral
#' Deconvolution and Identification System)
#' 
#' 
#' \code{parseELU} will typically be called by \code{\link{addAMDISPeaks}}, not
#' called directly.
#' 
#' Peaks that are detected within \code{rt.cut} are merged together.  This
#' avoids peaks which are essentially overlapping.
#' 
#' Fragments that are less than \code{min.pc} of the maximum intensity fragment
#' are discarded.
#' 
#' @param f ELU filename to read.
#' @param min.pc minimum percent of maximum intensity.
#' @param mz vector of mass-to-charge bins of raw data table.
#' @param rt.cut the difference in retention time, below which peaks are merged
#' together.
#' @param rtrange retention time range to parse peaks from, can speed up
#' parsing if only interested in a small region (must be \code{numeric} vector
#' of length 2)
#' @return
#' 
#' \code{list} with components \code{peaks} (table of spectra -- rows are
#' mass-to-charge and columns are the different detected peaks) and \code{tab}
#' (table of features for each detection), according to what is stored in the
#' ELU file.
#' @author Mark Robinson
#' @seealso \code{\link{addAMDISPeaks}}
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
#' eluFiles<-dir(gcmsPath,"ELU",full=TRUE)
#' 
#' # parse ELU file
#' eluList<-parseELU(eluFiles[1])
#'
#' @importFrom utils read.table
#' @importFrom stats cutree dist hclust median
#' @export parseELU
parseELU <- function(f, min.pc=0.01, mz=seq(50, 550), rt.cut=0.008, rtrange=NULL){
##    options(warn=-1)
    mostart <- function(u, key = "MO") {
        which(substr(u, 1, 2) == key)
    }
    v <- read.table(f, comment.char="", sep="\n", stringsAsFactors=FALSE)
    keep <- substr(v[, 1], 1, 5) == "NAME:"
    hdr <- v[keep, ]
    hdr <- strsplit(hdr, "\\|")
    rtis <- unlist(lapply(hdr, FUN=mostart, key="RT"))
    rts <- rep(NA, length(rtis))
    for(i in 1:length(rts)) rts[i] <- hdr[[i]][rtis[i]]
    rts <- as.numeric(gsub("RT", "", rts))
    starts <- which(substr(v[, 1], 1, 9) == "NUM PEAKS") + 1
    ends <- c(which(substr(v[, 1], 1, 5) == "NAME:")[-1] - 1, 
              length(v[, 1]))
    if(length(rtrange) == 2){
        w <- which(rts >= rtrange[1] & rts <= rtrange[2])
        rts <- rts[w]
        starts <- starts[w]
        ends <- ends[w]
    }else{
        w <- seq(1, length(starts))
    }
    peaks <- matrix(0, nrow=length(mz), ncol=length(starts))
    for(j in 1:length(starts)){## j <- 37 
        pks <- strsplit(paste(v[starts[j]:ends[j], 1], collapse=""), "\\)\\(")[[1]]
        pks[1] <- sub("\\(", "", pks[1])
        pks[length(pks)] <- sub("\\)", "", pks[length(pks)])
        aft <- lapply(pks, function(u, split=" ") .subset2(strsplit(u, split), 1)[2])
        aft <- unlist(lapply(aft, is.na))
        pks <- lapply(pks, function(u, split=" ") .subset2(strsplit(u, split), 1)[1])
        pks <- pks[aft]
        mzc <- as.numeric(unlist(lapply(pks, function(u, split=",") .subset2(strsplit(u, split), 1)[1])))
        idx <- which(mzc < range(mz)[1] | mzc > range(mz)[2])# is
                                        # there some m/z out of m/z
                                        # range of the experiment?
        if(length(idx) == 0){
            abn <- as.numeric(unlist(lapply(pks, function(u, split=",") .subset2(strsplit(u, split), 1)[2])))
        }else{
            mzc <- mzc[-idx] #
            abn <- as.numeric(unlist(lapply(pks, function(u, split=",") .subset2(strsplit(u, split), 1)[2])))
            abn <- abn[-idx] #
        }
        mzc <- mzc[abn > min.pc * max(abn)]
        abn <- abn[abn > min.pc * max(abn)]
        mt <- match(mzc, mz)
        mt <- mt[!is.na(mt)]
        ## cat(j, '\n')
        peaks[mt, j] <- abn        
    }
    v <- v[keep, ][w]
    v <- strsplit(v, "\\|")
    rts <- unlist(lapply(v, FUN=mostart, key="RT"))
    scs <- unlist(lapply(v, FUN=mostart, key="SC"))
    frs <- unlist(lapply(v, FUN=mostart, key="FR"))
    ras <- unlist(lapply(v, FUN=mostart, key="RA"))
    tab <- data.frame(matrix(NA, nrow=length(scs), ncol=5))
    colnames(tab) <- c("SC", "RT", "start", "end", "RA")
    for(i in 1:nrow(tab)){
        tab[i, 1] <- as.numeric(sub("SC", "", v[[i]][scs[i]]))
        tab[i, 2] <- as.numeric(sub("RT", "", v[[i]][rts[i]]))
        tab[i, 3:4] <- as.numeric(strsplit(sub("FR", "", v[[i]][frs[i]]), "-")[[1]])
        tab[i, 5] <- as.numeric(sub("RA", "", v[[i]][ras[i]]))
    }
    groups <- cutree(hclust(dist(tab[, 2])), h=rt.cut)
    newpeaks <- matrix(0, nrow=length(mz), ncol=max(groups))
    newtab <- data.frame(matrix(0, nrow=max(groups), ncol=5))
    colnames(newtab) <- colnames(tab)
    for(i in 1:nrow(newtab)){
        newtab[i, ] <- apply(tab[groups == i, ], 2, median)
        newtab[i, ncol(newtab)] <- sum(tab[groups == i, ncol(newtab)])
        newpeaks[, i] <- apply(matrix(peaks[, groups == i], nrow=length(mz)), 
            1, sum)
    }
    list(peaks=newpeaks, tab=newtab)
}





#' Parser for ChromaTOF files
#' 
#' Reads ASCII ChromaTOF-format files from AMDIS (Automated Mass Spectral
#' Deconvolution and Identification System)
#' 
#' 
#' \code{parseChromaTOF} will typically be called by
#' \code{\link{addChromaTOFPeaks}}, not called directly.
#' 
#' Peaks that are detected within \code{rt.cut} are merged together.  This
#' avoids peaks which are essentially overlapping.
#' 
#' Fragments that are less than \code{min.pc} of the maximum intensity fragment
#' are discarded.
#' 
#' @param fn ChromaTOF filename to read.
#' @param min.pc minimum percent of maximum intensity.
#' @param mz vector of mass-to-charge bins of raw data table.
#' @param rt.cut the difference in retention time, below which peaks are merged
#' together.
#' @param rtrange retention time range to parse peaks from, can speed up
#' parsing if only interested in a small region (must be \code{numeric} vector
#' of length 2)
#' @param skip number of rows to skip at beginning of the ChromaTOF
#' @param rtDivide multiplier to divide the retention times by (default: 60)
#' @return
#' 
#' \code{list} with components \code{peaks} (table of spectra -- rows are
#' mass-to-charge and columns are the different detected peaks) and \code{tab}
#' (table of features for each detection), according to what is stored in the
#' ChromaTOF file.
#' @author Mark Robinson
#' @seealso \code{\link{addAMDISPeaks}}
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
#' tofFiles<-dir(gcmsPath,"tof",full=TRUE)
#' 
#' # parse ChromaTOF file
#' cTofList<-parseChromaTOF(tofFiles[1])
#'
#' @importFrom utils read.table
#' @export parseChromaTOF
parseChromaTOF<-function(fn,min.pc=.01,mz=seq(85,500),rt.cut=.008,rtrange=NULL,skip=1, rtDivide=60) {
  f<-read.table(fn,sep="\t",quote="",comment.char="",skip=skip,header=TRUE,stringsAsFactors=FALSE)
  pk<-matrix(0,nrow=length(mz),ncol=nrow(f))
  for(i in 1:nrow(f)) {
    sp<-sapply(strsplit(f$Spectra[i]," ")[[1]],strsplit,split=":")
	pmz<-as.numeric(sapply(sp,.subset,1))
	int<-as.numeric(sapply(sp,.subset,2))
	m<-match(pmz,mz)
	pk[m,i]<-int
  }
  if( is.null(rtrange) )
    w <- seq_len( nrow(f) )
  else
    w <- which( (f$R.T/rtDivide) >= rtrange[1] & (f$R.T/rtDivide) <= rtrange[2] )
  list(tab=data.frame(rt=f$R.T[w],ht=f$Height[w]),peaks=pk[,w])
}
