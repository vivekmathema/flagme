
parseELU<-function(f,min.pc=.01,mz=seq(50,550),rt.cut=.008,rtrange=NULL) {
  mostart<-function(u,key="MO") { which(substr(u,1,2)==key) }
  v<-read.table(f,comment.char="",sep="\n",stringsAs=F)
  keep<-substr(v[,1],1,5)=="NAME:"
  hdr<-v[keep,]  # header lines
  hdr<-strsplit(hdr,"\\|")
  rtis<-unlist(lapply(hdr,FUN=mostart,key="RT"))
  rts<-rep(NA,length(rtis))
  for(i in 1:length(rts))
    rts[i]<-hdr[[i]][rtis[i]]
  rts<-as.numeric(gsub("RT","",rts))
  
  starts<-which(substr(v[,1],1,9)=="NUM PEAKS")+1
  ends<-c(which(substr(v[,1],1,5)=="NAME:")[-1]-1,length(v[,1]))
  if(length(rtrange)==2) {
    w<-which(rts>=rtrange[1] & rts<=rtrange[2])
    rts<-rts[w]
	starts<-starts[w]
	ends<-ends[w]
  } else {
    w<-seq(1,length(starts))
  }
  peaks<-matrix(0,nrow=length(mz),ncol=length(starts))
  for(j in 1:length(starts)) {
    pks<-strsplit(paste(v[starts[j]:ends[j],1],collapse=""),"\\)\\(")[[1]]
    pks[1]<-sub("\\(","",pks[1])
    pks[length(pks)]<-sub("\\)","",pks[length(pks)])
    aft<-lapply(pks,FUN=function(u,split=" ") .subset2(strsplit(u,split),1)[2])
	aft<-unlist(lapply(aft,is.na))
    pks<-lapply(pks,FUN=function(u,split=" ") .subset2(strsplit(u,split),1)[1])
	pks<-pks[aft]
    mzc<-as.numeric(unlist(lapply(pks,FUN=function(u,split=",") .subset2(strsplit(u,split),1)[1])))
    abn<-as.numeric(unlist(lapply(pks,FUN=function(u,split=",") .subset2(strsplit(u,split),1)[2])))
    mzc<-mzc[abn>min.pc*max(abn)]
    abn<-abn[abn>min.pc*max(abn)]
    peaks[match(mzc,mz),j]<-abn
  }
  v<-v[keep,][w]
  v<-strsplit(v,"\\|")
  rts<-unlist(lapply(v,FUN=mostart,key="RT"))
  scs<-unlist(lapply(v,FUN=mostart,key="SC"))
  frs<-unlist(lapply(v,FUN=mostart,key="FR"))
  ras<-unlist(lapply(v,FUN=mostart,key="RA"))
  tab<-data.frame(matrix(NA,nrow=length(scs),ncol=5))
  colnames(tab)<-c("SC","RT","start","end","RA")
  for(i in 1:nrow(tab)) {
    tab[i,1]<-as.numeric(sub("SC","",v[[i]][scs[i]]))
    tab[i,2]<-as.numeric(sub("RT","",v[[i]][rts[i]]))
	tab[i,3:4]<-as.numeric(strsplit(sub("FR","",v[[i]][frs[i]]),"-")[[1]])
    tab[i,5]<-as.numeric(sub("RA","",v[[i]][ras[i]]))
  }
  keep<-which(peaks[which(mz==73),]>0 | peaks[which(mz==147),]>0)  # remove peaks that have no evidence of TMS (73,147)
  tab<-tab[keep,]
  peaks<-peaks[,keep]
  groups<-cutree(hclust(dist(tab[,2])),h=rt.cut)  # merge peaks that have VERY similiar retention times
  newpeaks<-matrix(0,nrow=length(mz),ncol=max(groups))
  newtab<-data.frame(matrix(0,nrow=max(groups),ncol=5)); colnames(newtab)<-colnames(tab)
  for(i in 1:nrow(newtab)) {
    newtab[i,]<-apply(tab[groups==i,],2,median)
	newtab[i,ncol(newtab)]<-sum(tab[groups==i,ncol(newtab)])
    newpeaks[,i]<-apply(matrix(peaks[,groups==i],nrow=length(mz)),1,sum)
  }
  list(peaks=newpeaks,tab=newtab)
}

parseChromaTOF<-function(fn,min.pc=.01,mz=seq(85,500),rt.cut=.008,rtrange=NULL,skip=1) {
  f<-read.table(fn,sep="\t",quote="",comment.char="",skip=skip,header=TRUE,stringsAsFactors=FALSE)
  pk<-matrix(0,nr=length(mz),nc=nrow(f))
  for(i in 1:nrow(f)) {
    sp<-sapply(strsplit(f$Spectra[i]," ")[[1]],strsplit,split=":")
	pmz<-as.numeric(sapply(sp,.subset,1))
	int<-as.numeric(sapply(sp,.subset,2))
	m<-match(pmz,mz)
	pk[m,i]<-int
  }
  list(tab=data.frame(rt=f$R.T,ht=f$Height),peaks=pk)
}
