# summarizing results of GLM models and bootstraps on thresholds per group of species


summary.glmtest <- function(M = glmSR.overall,data=species, group="ALIEN", 
                            type =c("glm", "boot", "overall.boot"),threshold= "th.CI") {
  
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
    thresh.exist = matrix(!is.na(S$thresh[, threshold]),nrow =dim(S$dif)[1], ncol =dim(S$dif)[2], byrow=F)
    above.thr = as.data.frame(sapply(1:dim(S$dif)[2], function(x) (S$thresh[, threshold]-1)<=x), row.names=rownames(S$dif), names=names(S$dif))
    thr.is = as.data.frame(sapply(1:dim(S$dif)[2], function(x) (S$thresh[, threshold]-1)==x), row.names=rownames(S$dif), names=names(S$dif))
    
    # number of potentieally signif species occurrence per class
    n.sp <- apply(!is.na(S$dif), 2, sum, na.rm=T)
    # number of species which have a threshold per class
    n.target    <-  apply((!is.na(S$dif) & !is.na(thresh.exist)), 2, sum, na.rm=T)
    
    
    # number of significant negative effects per class
    n.negative  <-  apply(S$z<0, 2, sum, na.rm=T)
  
    # number of significant negative effects per class
    n.negative.target<-apply( (S$z<0 & !is.na(S$dif) & above.thr), 2, sum, na.rm=T)
         
    
    # number of significant positive effects per class
    n.positive<-sapply(names(S$dif), FUN= function(x) {
      length(which(S$z[,x]>=0))
    })
    
    # number of significant negative impavts per class for species with a threshld
    n.impact<-sapply(names(S$dif), FUN= function(x) {
      length(which(S$z[,x]<0 & !is.na(S$thresh[ , threshold])))
    })
    
    # proportion of signifi negative effects per class
    p.impact<-n.impact/ n.sp
    
    # proportion of signifi negative effects per class
    p.negative<-n.negative/ n.sp
    
    # number of times each class is the threshold
#     freq.thr <-table(as.factor(S$thresh[, threshold]))[match(2:6,names(table(as.factor(S$thresh[, threshold]))))]
    freq.thr <- apply(thr.is, 2, sum, na.rm=T)
    names(freq.thr)<-names(n.sp)
    freq.thr[is.na(freq.thr)]<-0

    prop.thr <- freq.thr / n.sp
    if(sum(freq.thr, na.rm=T)!=0) perc.thr <- freq.thr / sum(freq.thr, na.rm=T)
    if(sum(freq.thr, na.rm=T)==0) perc.thr <- freq.thr 
    
    out<-rbind(out, data.frame(group=names[j],class = 2:6, nb.sp=n.sp, n.target = n.target,n.negative.target = n.negative.target, freq.impact=n.impact, freq.negative =n.negative,freq.positive =n.positive,
                               prop.impact=as.array(p.impact),freq.thr=freq.thr, prop.thr=prop.thr, perc.thr=perc.thr))
    
    #summary table
    n.sp <- length(rownames(S$thresh))
    freq.thr <- length(which(!is.na(S$thresh[, threshold])))
    
    sum.table <-rbind(sum.table,  data.frame(group=names[j], nb.sp=n.sp,  freq.thr=freq.thr))
    
  }
  return(list(class.summary =out, overall.summary = sum.table))
}
