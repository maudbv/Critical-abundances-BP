
pairedclass.test=function(db=databp, var='SR',domclass='domlevels',zeros=F,env=envplot, min.occur=15,  min.class=2, alpha=0.05) {
 
  #  function testing all relationships between dominance leveles (class) and a variable of the community (e.g. SR)

  db$Sp.occurence=table(as.character(db$SpeciesCode))[match(as.character(db$SpeciesCode),unique(as.character(db$SpeciesCode)))]
  tmp=db[db$Sp.occurence>min.occur,] 
  
  # initiate results
  ts= data.frame(c2=NA, c3=NA, c4=NA, c5=NA, c6=NA)
  dfs=data.frame(c2=NA, c3=NA, c4=NA, c5=NA, c6=NA)
  ps=data.frame(c2=NA, c3=NA, c4=NA, c5=NA, c6=NA)
  kruskal=data.frame(k=NA, df=NA, p.val=NA)
  anov=data.frame(F.val=NA, df=NA, p.val=NA)
  diff=data.frame(c2=NA, c3=NA, c4=NA, c5=NA, c6=NA)
  
  # looping on species
  for (i in 1:length(unique(tmp$SpeciesCode)) ) {
    
    sp=unique(tmp$SpeciesCode)[i]
    x=tmp[as.character(tmp$SpeciesCode)==sp,c(domclass,var,'PlotName')]
    names(x)=c('abun','var','PlotName')
    z= cbind(abun=NA,env[which(!env$PLOTID %in% x$PlotName),c(var,"PLOTID")])
    x=apply(x,2,as.numeric)
    z=apply(z,2,as.numeric)
    x=as.data.frame(rbind(x,z))

    x$abun=7-as.numeric(as.character(x$abun))
    if (zeros) x$abun[is.na(x$abun)]=0
    
    abun=sort(as.numeric(as.character(na.omit(unique(x$abun)))))
  
    print(sp)
    # if the species is not represented in at least 'min.abun' different levels of abundance :
#     if(length(unique(x$abun)) < min.class) {
#       r=c(unique(tmp$Sp.occurence[i]),nabun,
#       max(x$abun, na.rm=T), min(x$abun[x$abun>0], na.rm=T),
#       NA,NA, NA,NA,NA)
#       }
#       print(length(unique(x$abun)))
    
    # if the species is represented in sufficient different levels of abundance :
    if(length(unique(x$abun))>=min.class) {
      
        # kruskal-wallis test across all classes
        k=kruskal.test(x$var,x$abun)  
        kruskal[i,] =  c(k$stat, k$param, k$p.value)
        
        # one way anova test across all classes
        a=anova(lm(var~ as.factor(abun), data=x))
        anov[i,] =  c(a$F[1], a$Df[1]+a$Df[2], a$P[1])
      
        # contrast min abund to all 5 others
        m=min(abun)
        j=m
        low=x[which(x$abun==m),"var"]
        while (j<max(abun)) {
          high=x[which(x$abun==j+1),"var"]
           if (length(high)>5) { 
             r=t.test(low,high)
             ts[i,j] =  r$stat
             dfs[i,j] = r$param
             ps[i,j] =  r$p.value
             diff[i,j] =  r$est[2] -r$est[1]
           }
            j=j+1
          
            print(j)
         }
      }
    }
  rownames(kruskal) <- rownames(anov) <-rownames(ts) <- rownames(dfs) <-rownames(diff) <-rownames(ps) <-unique(tmp$SpeciesCode)
  
#### thresholds
    th=data.frame(t( sapply(rownames(ps), FUN=function(y)  {
      o=rep(NA, length(names(ps)))
      tmp=which(!is.na(ps[y,]) & ps[y,]<=0.05 & ts[y,]>0)
      tmp2=which(!is.na(ts[y,]) & ts[y,]>0)
      if(length(tmp)>0 ) { if (all(min(tmp):max(tmp)%in% tmp2))   o[min(tmp)]=1}
      return(o)
    })))
  
    names(th)=names(ps)
    threshold=data.frame(threshold=apply(th,1, FUN=function(x) {
      t=NA
      if (length(which(x==1))>0) t=which(x==1)
      return(t)
      }))
browser()
  
 ### % of plots over threshold
 nplots= t(sapply(rownames(threshold), function(i) {
   n=as.vector(rep(NA,6))
   names(n) = 1:6
   o=table(databp[databp$SpeciesCode == i, "domlevels"])
   n[names(o)]  = o 
   return(n)
 }))
  
threshold$nplots.impact = sapply(rownames(threshold), function(i) {
   length(databp[databp$SpeciesCode == i, "domlevels"] >threshold[i,"threshold"])
 })
  
threshold$pplots.impact = sapply(rownames(threshold), function(i) {
    length(which(databp[databp$SpeciesCode == i, "domlevels"] >threshold[i,"threshold"]))/length(databp[databp$SpeciesCode == i, "domlevels"])
  })

  
  return(list(kruskal.test=kruskal,anova=anov, t=ts, df=dfs, p.val=ps, threshold=threshold, diff=diff))
  
  
}


plot.pc=function(pc=pclassSR) {
  par(mfrow=c(1,3))
  # number of signifi negative effects per class
  barplot(sapply(names(pclassSR$p.val), FUN= function(x) {
    length(which(pclassSR$p.val[,x]<=0.05 &pc$t[,x]>0))
  }))
  title(ylab="number of negative effects")
  
  # proportion of signifi negative effects per class
  barplot(sapply(names(pclassSR$p.val), FUN= function(x) {
    length(which(pclassSR$p.val[,x]<=0.05 &pc$t[,x]>0))/length(which(!is.na(x) &pc$t[,x]>0))
  }))
  title(ylab="% of negative effects")
  
  # number of times each class is the minimum abundance for effects
  barplot((f <- function(X=pclassSR) {
    tmp=data.frame(t( sapply(rownames(pclassSR$p.val), FUN=function(y) {  
      o=rep(NA, length(names(pclassSR$p.val)))
      tmp=which(!is.na(pclassSR$p.val[y,]) &pc$p.val[y,]<=0.05 &pc$t[y,]>0)
      if(length(tmp)>0 )  o[min(tmp)]=1
      return(o)
    })))
    names(tmp)=names(pclassSR$p.val)
    p=apply(tmp, 2, sum, na.rm=T)
  })())
  title(ylab="Freq. of first effect")
}

  