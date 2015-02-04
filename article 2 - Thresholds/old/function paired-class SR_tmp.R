
pairedclass.test=function(db=databp, var='SR',domclass='abun',zeros=F, verbose=F, test="w",
                          env=envplot, min.occur=15,  min.class=2, alpha=0.05) {
  require(coin)
  #  function testing difference in a community metric (e.g. SR) between low abundances and increasing abundance classes
  
  # calculate species prevalence
  db$Sp.occurence=table(as.character(db$SpeciesCode))[match(as.character(db$SpeciesCode),unique(as.character(db$SpeciesCode)))]
  tmp=db[db$Sp.occurence>min.occur,] 
  
  # initiate results
  ts= data.frame(c2=NA, c3=NA, c4=NA, c5=NA, c6=NA)
  dfs=data.frame(c2=NA, c3=NA, c4=NA, c5=NA, c6=NA)
  ps=data.frame(c2=NA, c3=NA, c4=NA, c5=NA, c6=NA)
  kruskal=data.frame(k=NA, df=NA, p.val=NA)
  anov=data.frame(F.val=NA, df=NA, p.val=NA)
  diff=data.frame(c2=NA, c3=NA, c4=NA, c5=NA, c6=NA)
  n.obs = data.frame(c2=NA, c3=NA, c4=NA, c5=NA, c6=NA)
  
  # looping on species
  for (i in 1:length(unique(tmp$SpeciesCode)) ) {
    
    sp=unique(tmp$SpeciesCode)[i]
    x=tmp[as.character(tmp$SpeciesCode)==sp,c("abun",var,'PlotName')]
    names(x)=c('abun','var','PlotName')
    z= cbind(abun=NA,env[which(!env$PLOTID %in% x$PlotName),c(var,"PLOTID")])
    x=apply(x,2,as.numeric)
    z=apply(z,2,as.numeric)
    x=as.data.frame(rbind(x,z))
    abun=sort(as.numeric(as.character(na.omit(unique(x$abun)))))
    
    if (verbose) print(sp)
    
    if(length(na.omit(unique(x$abun)))>=min.class) {
      
      # kruskal-wallis test across all classes
      k=kruskal.test(x$var,x$abun)  
      kruskal[i,] =  c(k$stat, k$param, k$p.value)
      
      # one way anova test across all classes
      a=anova(lm(var~ as.factor(abun), data=x))
      anov[i,] =  c(a$F[1], a$Df[1]+a$Df[2], a$P[1])
      
      # contrast min abund to all others
      m=min(abun)
      j=m
      low=x[which(x$abun==m),"var"]
      while (j<max(abun)) {
        high=x[which(x$abun==j+1),"var"]
        n.obs[i,j] =  length(high)

        diff[i,j] =  mean(high, na.rm=T) - mean (low, na.rm=T)
        if (length(high)>3) {    ### more than 3 observations in the abundance class

          if (test=="t") {
            r=t.test(low,high)
            ts[i,j]  =  r$stat
            dfs[i,j] =  r$param
            ps[i,j]  =  r$p.value
          }
          
          if (test=="w"){
            d=data.frame(sr=c(high,low), group=c(rep(1, length(high)), rep(2, length(low))))
            r=wilcox_test(sr ~ as.factor(group), data=d)
            ps[i,j] =  pvalue(r)
            ts[i,j]  =  NA
            dfs[i,j] =  NA
          }
        }
        j=j+1
       if (verbose)  print(j)
      }
    }
  }

  rownames(kruskal) <- rownames(anov) <-rownames(ts) <- rownames(dfs) <-rownames(diff) <-rownames(ps) <-rownames(n.obs) <-unique(tmp$SpeciesCode)

  #### Detecting thresholds of negative impact
  
  th=data.frame(t( sapply(rownames(ps), FUN=function(y)  {
    # presence or absence of abundance threshold for a species:
    o=rep(NA, length(names(ps)))
    max=NA
#     if (length(which(!is.na(ts[y,])))>0) max=max(which(!is.na(ts[y,])), na.rm=T)
#     tmp=which(!is.na(ps[y,]) & ps[y,]<=0.05 & diff[y,]<0) ## alpha = 5%
#     tmp2=which(!is.na(ts[y,]) & diff[y,]<0)  ## subsequent classes have negative differences, ans sufficient statistical power for the test
#     if(length(tmp)>0 ) { if (all(min(tmp):max %in% tmp2))   o[min(tmp)]=1}
    
    if (length(which(!is.na(ps[y,])))>0) max=max(which(!is.na(diff[y,])), na.rm=T)
    tmp=which(!is.na(ps[y,]) & ps[y,]<=0.05 & diff[y,]<0) ## alpha = 5%
    tmp2=which(diff[y,]<0)  ## subsequent classes have negative differences, no matter significance or sample size??
    if(length(tmp)>0 ) { if (all(min(tmp):max %in% tmp2))   o[min(tmp)]=1}
    return(o)
  })))

  
  names(th)=names(ps)
  threshold=data.frame(threshold=apply(th,1, FUN=function(x) {
    t=NA
    if (length(which(x==1))>0) t=which(x==1)+1
    return(t)
  }))
  
  nplots= t(sapply(rownames(threshold), function(i) {
    n=as.vector(rep(NA,6))
    names(n) = 1:6
    o=table(db[db$SpeciesCode == i, "abun"])
    n[names(o)]  = o 
    return(n)
  }))
  
  ### % of plots over threshold
  threshold$nplots.impact = sapply(rownames(threshold), function(i) {
    x=NA  
    if (!is.na(threshold[i,"threshold"])) x= sum(db[db$SpeciesCode == i, "abun"] >= threshold[i,"threshold"], na.rm=T)
    return(x)
  })
  
  threshold$pplots.impact = sapply(rownames(threshold), function(i) {
    x=NA  
    if (!is.na(threshold[i,"threshold"])) {
      x=sum(db[db$SpeciesCode == i, "abun"] >= threshold[i,"threshold"], na.rm=T)/length(na.omit(db[db$SpeciesCode == i, "abun"]))
    }
    return(x)
  })
  
  #### impact sizes
  threshold$mean.sig.diff = sapply(rownames(threshold), function(i) {
    x=NA  
    if (!is.na(threshold[i,"threshold"])) x= mean(as.numeric(diff[i, which(1:5 >= threshold[i,"threshold"] -1)]), na.rm=T)
    return(x)
  })
  
  threshold$wtd.mean.sig.diff = sapply(rownames(threshold), function(i) {
    x=NA  
    if (!is.na(threshold[i,"threshold"])) {
      d=as.numeric(diff[i, which(1:5 >=threshold[i,"threshold"] - 1)])
      nb=as.numeric(n.obs[i, which(1:5 >=threshold[i,"threshold"] - 1)])
      x=sum(d*nb, na.rm=T)/sum(nb, na.rm=T)
    }
    return(x)
  })
  
  threshold$max.sig.diff = sapply(rownames(threshold), function(i) {
    x=NA  
    if (!is.na(threshold[i,"threshold"])) x= min(as.numeric(diff[i, which(1:5 >=threshold[i,"threshold"]-1)]), na.rm=T)
    return(x)
  })
  
  threshold$thr.diff = sapply(rownames(threshold), function(i) {
    x=NA  
    if (!is.na(threshold[i,"threshold"])) x= as.numeric(diff[i, threshold[i,"threshold"]-1])
    return(x)
  })
  
  threshold$prevalence = sapply(rownames(threshold), function(i) {length( db[db$SpeciesCode == i, "abun"] ) })
  
  return(list(kruskal.test=kruskal,anova=anov, t=ts, df=dfs, p.val=ps, diff=diff,n.obs=n.obs, threshold=threshold ))
}


plot.sp=function(sp="ULEEUR", pc=pclassSR, var="SR", db=databp)
  {  
      Z= db[db$SpeciesCode==sp,]
      boxplot(as.formula(paste(var, " ~ abun")), data=Z, xlim=c(0,6), outline=F, varwidth=T)
      mtext(side=3, line=-1.1, at=2:6,text=sapply(pc$p.val[sp,], FUN=p2star), cex=0.8)
      mtext(side=3, text=paste("K:",p2star(pc$kruskal.test[sp,"p.val"])), adj=1)
      title(main=sp, cex.main=0.8)
}

plot.pc=function(pc=pclassSR, var="SR", db=databp, sel.criteria = c("kruskal", "anova", "spearman", "none"),
                 spear = effect.summary) {
  par(mfrow=c(5,5), mar=c(2,2,2,2))
  a=0
  for (i in rownames(pc$p.val)){
    # only species with significant kruskal walis test :
    if (sel.criteria == "kruskal") sel <- pc$kruskal.test[i,"p.val"] <=0.05
    if (sel.criteria == "anova") sel <- pc$anova[i,"p.val"] <=0.05
    if (sel.criteria == "anova") sel <- spear[i,"S.pval"] <=0.05
    if (sel.criteria == "none") sel <- TRUE
    
    if (sel==T){
      plot.sp(sp= i, pc=pc, var=var, db=db)
      a=a+1
    }
  }
}

### Summarize results :
summary.pc=function(pc=pclassSR, data=species, group="ALIEN", graph=T) {
    data=data[rownames(pc$p.val),]
    if (is.null(group)) {
    sub=pc
    G=1
  }
  
  if (!is.null(group)) { 
    sub=list()
    i=1
    for (g in sort(unique(data[,group]))) {
      sub[[i]]=lapply(pc, FUN= function(x) x[rownames(x) %in% rownames(data)[data[,group]==g],] )
      i=i+1
      G=length(unique(data[,group]))
      names=paste(group, ":",sort(unique(data[,group])), sep="")
    }
  }                                 
   
   par(mfrow=c(G,2), oma=c(4,6,6,2), mar=c(2,2,1,1))
  
  out=NULL   
  for (j in 1:G) {
    S=sub[[j]]
    
    # number of potentieally signif species occurrence per class
    n.sp=sapply(names(S$p.val), FUN= function(x) {
      length(which(!is.na(S$p.val[,x])))
    })
    
    # number of significant negative effects per class
    n.impact=sapply(names(S$p.val), FUN= function(x) {
      length(which(S$p.val[,x]<=0.05 & S$diff[,x]<0))
    })
    
    # number of significant negative effects per class
    n.positives=sapply(names(S$p.val), FUN= function(x) {
      length(which(S$p.val[,x]<=0.05 & S$diff[,x]>0))
    })
    
    # proportion of signifi negative effects per class
    p.impact=n.impact / colSums(!is.na(S$p.val))
    # number of times each class is the threshold
    freq.thr= (f <- function(X=SlassSR) {
      tmp=table(S$threshold$threshold)
      y=as.vector(rep(0,5)) 
      names(y)=names(S$p.val)
      y[as.numeric(names(tmp))-1]=tmp
      return(y)
    })()

    # proportion of times a class is a threshold
    prop.thr = freq.thr / colSums(!is.na(S$p))
    perc.thr = freq.thr / sum(freq.thr, na.rm=T)
    
    #graphical representation
 if (graph) {
    barplot( p.impact)
    if (j==1) mtext(side=3, text="% negative\neffects",adj=0.5,line=1, cex=0.8)
    
    barplot(prop.thr, col="black" )
    if (j==1) mtext(side=3, text="% threshold\neffects", adj=0.5,line=1, cex=0.8)
    
    out=rbind(out, data.frame(group=names [j],nb.sp=n.sp, freq.impact=n.impact, freq.positives =n.positives,
                              prop.impact=p.impact,freq.thr=freq.thr, prop.thr=prop.thr, perc.thr=perc.thr))
  }
  

  mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=3, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
  mtext(text="Abundance class", side=1, outer=T, line=2)
  }
  return(out)
}



########## MORE spline than density : to be improved

threshold.spline=function(data=effects.classSR,index="prop.thr", m= "natural",
                          tit="Total richness", ylab="Proportion of threshold per class") {
  y1 <- data[data$group=="ALIEN:0", index]
  y2 <- data[data$group=="ALIEN:1", index]
  f1 <- splinefun(2:6,y1, method=m)
  f2 <- splinefun(2:6,y2, method=m)
  x=seq(2,6, 0.01)
  
  plot(x,f1(x), type="l", ylim=c(0,max(c(y1,y2))+0.05), ann=F, las=1)
  lines(x,f2(x), col="red")
  points(2:6,y1, type="p")
  points(2:6,y2, type="p",col="red")
  
  mtext(side=1, "Abundance class", line= 2)
  mtext(side=2, text=ylab, line= 3.5)
  mtext(side=3, tit, line= 1)
  
  legend("topleft", legend=c("Alien targets", "Native targets"), col=c("red", "black"), lty="solid", bty="n")
}