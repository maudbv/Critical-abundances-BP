# summarizing results of GLM models and bootstraps on thresholds per group of species


summary.glmtest <- function(M = glmSR,data=species, group="ALIEN", type =c("glm", "boot", "overall.boot")) {
  
  ### select onlys species in the model
  data=data[rownames(M$dif),]
  
  if (type == "boot") M$thresh <- M$boot.thresh
  if (type == "overall.boot") M$thresh <- M$impact.spread
  
  ## divide species according to grouping factor (ie Alien vs Native)
  if (is.null(group)) {
    sub<-M
    G<-1
  }
  
  if (!is.null(group)) { 
    sub<-list()
    i<-1
    G<-length(unique(data[,group]))
    for (g in sort(unique(data[,group]))) {
      sub[[i]]<-lapply(M, FUN= function(x) x[rownames(x) %in% rownames(data)[data[,group]==g],] )
      i<-i+1
    }
    names<-paste(group, ":",sort(unique(data[,group])), sep="")
  }                                 
  
  out<-NULL   
  sum.table <- NULL 
  
  for (j in 1:G) {
    S<-sub[[j]]
    
    # number of potentieally signif species occurrence per class
    n.sp<-sapply(names(S$dif), FUN= function(x) {
      length(which(!is.na(S$dif[,x])))
    })
    
    # number of species which have a threshold per class
    n.target<-sapply(names(S$dif), FUN= function(x) {
      length(which(!is.na(S$dif[,x]) & !is.na(S$thresh$th)) )
    })
    
    
    # number of significant negative effects per class
    n.negative<-sapply(names(S$dif), FUN= function(x) {
      length(which(S$P[,x]<=0.05 & S$z[,x]<0))
    })
    
    
    # number of significant positive effects per class
    n.positive<-sapply(names(S$dif), FUN= function(x) {
      length(which(S$P[,x]<0.05 & S$z[,x]>=0))
    })
    
    # number of significant negative impavts per class for species with a threshld
    n.impact<-sapply(names(S$dif), FUN= function(x) {
      length(which(S$P[,x]<=0.05 & S$z[,x]<0 & !is.na(S$thresh$th)))
    })
    
    # proportion of signifi negative effects per class
    p.impact<-n.impact/ n.sp
    
    # proportion of signifi negative effects per class
    p.negative<-n.negative/ n.sp
    
    # number of times each class is the threshold
    freq.thr <-table(as.factor(S$thresh$th))[match(2:6,names(table(as.factor(S$thresh$th))))]
    names(freq.thr)<-names(n.sp)
    freq.thr[is.na(freq.thr)]<-0
    prop.thr <- freq.thr / n.sp
    if(sum(freq.thr, na.rm=T)!=0) perc.thr <- freq.thr / sum(freq.thr, na.rm=T)
    if(sum(freq.thr, na.rm=T)==0) perc.thr <- freq.thr 
    
    out<-rbind(out, data.frame(group=names[j],class = 2:6, nb.sp=n.sp, n.target = n.target, freq.impact=n.impact, freq.negative =n.negative,freq.positive =n.positive,
                               prop.impact=as.array(p.impact),freq.thr=freq.thr, prop.thr=prop.thr, perc.thr=perc.thr))
    
    #summary table
    n.sp <- length(rownames(S$thresh))
    freq.thr <- length(which(!is.na(S$thresh$th)))
    
    sum.table <-rbind(sum.table,  data.frame(group=names[j], nb.sp=n.sp,  freq.thr=freq.thr))
    
  }
  return(list(class.summary =out, overall.summary = sum.table))
}
