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
      sub[[i]] <- lapply(M,function(x) {   ## We subset each element of M by alien vs native
        if (!all(!rownames(x) %in% rownames(data)[data[,group]==g]) & !is.null(rownames(x))) 
        {
          y <-  x[rownames(x) %in% rownames(data)[data[,group]==g],] 
        }
        else y = NA
        return(y)
      })
            i<-i+1
    }
    group.names<-paste(group, ":",sort(unique(data[,group])), sep="")
  }

  out<-NULL
  sum.table <- NULL

  for (j in 1:G) {
    S<-sub[[j]]
    thresh.exist = matrix(!is.na(S$thresh[, threshold]),nrow =dim(S$dif)[1], ncol =dim(S$dif)[2], byrow=F)
    above.thr = as.data.frame(sapply(1:dim(S$dif)[2], function(x) (S$thresh[, threshold]-1)<=x), row.names=rownames(S$dif), names=names(S$dif))
    thr.is = as.data.frame(sapply(1:dim(S$dif)[2], function(x) (S$thresh[, threshold]-1)==x), row.names=rownames(S$dif), names=names(S$dif))

    testCI = as.data.frame(t(sapply(rownames(S$dif), FUN=function(l) (sign(M$CIlow[l,]) * sign(M$CIhi[l,]))==1)),
                             row.names=rownames(S$dif), names=names(S$dif))


    # number of potentieally signif species occurrence per class
    n.sp <- apply(!is.na(S$dif), 2, sum, na.rm=T)
    n.sp <- apply(!is.na(S$z), 2, sum, na.rm=T)

    # number of species which have a threshold per class
    n.target<-  apply((!is.na(S$z) & !is.na(thresh.exist)), 2, sum, na.rm=T)


    # number of significant positive effects per class
    n.positive<-apply( S$z>0 & !is.na(S$dif), 2, sum, na.rm=T)

    # number of negative effects per class
    n.negative  <-  apply(S$z<0 & !is.na(S$dif), 2, sum, na.rm=T)

    # number of negative effects significant according to bootstrap CI
    n.negative.sig<-apply( (S$z<0 & !is.na(S$dif) & testCI), 2, sum, na.rm=T)

    # number of negative effects significant according to bootstrap CI
    n.positive.sig<-apply( (S$z>0 & !is.na(S$dif) & testCI), 2, sum, na.rm=T)

    # number of negative effects above thresholds per class
    n.negative.above<-apply( (S$z<0 & !is.na(S$dif) & above.thr), 2, sum, na.rm=T)

       # proportion of signifi negative effects per class
    p.negative<-n.negative/ n.sp
    p.negative.sig <- n.negative.sig/ n.sp
    # number of times each class is the threshold
#     freq.thr <-table(as.factor(S$thresh[, threshold]))[match(2:6,names(table(as.factor(S$thresh[, threshold]))))]
    freq.thr <- apply(thr.is, 2, sum, na.rm=T)
    names(freq.thr)<-names(n.sp)
    freq.thr[is.na(freq.thr)]<-0

    prop.thr <- freq.thr / n.sp
    if(sum(freq.thr, na.rm=T)!=0) perc.thr <- freq.thr / sum(freq.thr, na.rm=T)
    if(sum(freq.thr, na.rm=T)==0) perc.thr <- freq.thr

    out<-rbind(out, data.frame(group=group.names[j],class = 2:6, nb.sp=n.sp, n.target = n.target,
                               freq.negative = n.negative, freq.positive =n.positive,
                               freq.negative.sig = n.negative.sig, freq.positive.sig =n.positive.sig,
                               freq.negative.above = n.negative.above,
                               prop.negative = p.negative,prop.negative.sig = p.negative.sig,
                               freq.thr=freq.thr, prop.thr=prop.thr, perc.thr=perc.thr))

    #summary table
    n.sp <- length(rownames(S$thresh))
    freq.thr <- length(which(!is.na(S$thresh[, threshold])))

    sum.table <-rbind(sum.table,  data.frame(group=group.names[j], nb.sp=n.sp,  freq.thr=freq.thr))

  }
  return(list(class.summary =out, overall.summary = sum.table))
}
