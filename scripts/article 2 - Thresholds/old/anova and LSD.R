### detect threshold with global ANOVa or spearman + then post hoc LSD test
LSDclass.test=function(db=databp, var='SR',domclass='abun',zeros=F, verbose=F, test="w",
                       env=envplot, min.occur=10,  min.class=3, alpha=0.05) {
  require(coin)

## sp spanning the minimum number of abundance classes
a <- unique(names(which(rowSums(table(db$SpeciesCode, db$abun)>0)>=min.class)))  

## species with minimum occurrences
# calculate species prevalence
db$Sp.occurence <- table(as.character(db$SpeciesCode))[match(as.character(db$SpeciesCode),unique(as.character(db$SpeciesCode)))]
b <- unique(names(which(db$Sp.occurence>min.occur)))

# select only species in a and b groups :
db.modif <- db[which( (db$SpeciesCode %in% a) & (db$SpeciesCode %in% b)),] 

# list of species name to be included in analysis :
sp.names <- unique(db.modif$SpeciesCode)


  # initiate results
  diff <- leastdiff <- n.obs <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 5,
                                                       dimnames=list(sp.names,c("c2", "c3", "c4", "c5", "c6"))))
  kruskal=data.frame(matrix(NA, nrow=length(sp.names), ncol= 3,
                            dimnames=list(sp.names, c("k" ,"df", "p.val"))))

  spear=data.frame(matrix(NA, nrow=length(sp.names), ncol= 3,
                          dimnames=list(sp.names, c("k" ,"df", "p.val"))))
  anov <- data.frame(matrix(NA, nrow=length(sp.names), ncol= 4,
                             dimnames=list(sp.names, c("F" ,"df", "MSerr", "p.val"))))

  # looping on species
  for (i in 1:length(sp.names) ) {
   
    sp=sp.names[i]  # select species name
    x=db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",var,'PlotName')] # select occurrences of the species
    names(x)=c('abun','var','PlotName')
#     z= cbind(abun=NA,env[which(!env$PLOTID %in% x$PlotName),c(var,"PLOTID")])
#     x=apply(x,2,as.numeric)
#     z=apply(z,2,as.numeric)
#     x=as.data.frame(rbind(x,z))
    
abun=sort(as.numeric(as.character(na.omit(unique(x$abun))))) ## list of abundance classes for species i
n.obs[i,]=sapply(2:6, FUN=function(j) length(x[which(x$abun==j),"var"]))

   
      
      #calculate difference in class mean SR
      diff[i,] =  sapply(2:6, FUN=function(j) mean(x[which(x$abun==j),"var"], na.rm=T) - mean (x[which(x$abun==min(abun)),"var"], na.rm=T))
      
      # kruskal-wallis test across all classes
      k=kruskal.test(x$var,x$abun)  
      kruskal[i,] =  c(k$stat, k$param, k$p.value)
      
      # spearmnan test across all classes
      s =cor.test(x$var,x$abun, method="spearman")  
      spear[i,] =  c(s$estimate, NA, s$p.value)
      
      # one way anova test across all classes
      a=aov(var~ as.factor(abun), data=x)
      anov[i,] =  c(anova(a)[1,4], anova(a)[2,1],anova(a)[2,3], anova(a)[1,5])
      
      if(anov[i,"p.val"]<=0.05) {
        # LSD post hoc test for all classes
#         l=LSD.test(a, "as.factor(abun)")
        tc <- qt(0.975, df= anov[i,"df"])
        leastdiff [i,] =sapply(1 : 5, FUN=function(j) {
          least= tc * sqrt(anov[i,"MSerr"] *(1/n.obs[i,1] + 1/n.obs[i,j]))
          return(least)
        })
      }

    }


rownames(kruskal) <- rownames(anov) <-rownames(spear) <-rownames(anov) <-rownames(leastdiff)<-rownames(diff) <-rownames(n.obs)<-sp.names 

#  return(list(n.obs=n.obs, anova=anov, kruskal= kruskal, diff=diff,leastdiff=leastdiff ))


##### Detecting thresholds of negative impact
th=data.frame(t( sapply(rownames(diff), FUN=function(y)  {
 th=NA # initiate : no threshold
 if(sum(diff[y,] +leastdiff[y,]<0, na.rm=T)>0) {  # check there is at least one difference to test
 m = min(which(diff[y,] +leastdiff[y,]<0), na.rm=T )   ### detect first significant difference
if (all(diff[y, min(length(names(diff)),(m+1)):length(names(diff))] <0, na.rm=T)) {## check all subsequent classes have negative differences
  th=m+1  # the threshold abundance class is m+1
}
}
return(th)
})))


### summarizing threshold information

threshold=data.frame(threshold=t(th))
### % of plots over threshold

nplots= t(sapply(rownames(threshold), function(i) {
  n=as.vector(rep(NA,6))
  names(n) = 1:6
  o=table(db[db$SpeciesCode == i, "abun"])
  n[names(o)]  = o 
  return(n)
}))

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
  if (!is.na(threshold[i,"threshold"])) x= min(as.numeric(diff[i, which(1:5 >=threshold[i,"threshold"]-1)]), na.rm=T) # maximum diff is in fact the minimum because negative values
  return(x)
})

threshold$thr.diff = sapply(rownames(threshold), function(i) {
  x=NA  
  if (!is.na(threshold[i,"threshold"])) x= as.numeric(diff[i, threshold[i,"threshold"]-1])
  return(x)
})

threshold$prevalence = sapply(rownames(threshold), function(i) {length( db[db$SpeciesCode == i, "abun"] ) })

return(list(kruskal.test=kruskal,anova=anov, spearman=spear, diff=diff,leastdiff=leastdiff, n.obs=n.obs, threshold=threshold ))

}


plot.sp.lsd=function(sp="ULEEUR", lsd=LSDclassSR, var="SR", db=databp)
{  
  Z= db[db$SpeciesCode==sp,]
  
  col=rep("white",6)
  col[ as.numeric( lsd$threshold[sp,"threshold"])]="red"
  
  boxplot(as.formula(paste(var, " ~ abun")), data=Z, xlim=c(0,6), outline=F, col=col,varwidth=T)
  
  
  sig = c("ns","*")[(lsd$diff[sp,] + lsd$leastdiff[sp,]<0) + 1]
  sig[is.na(sig)] = ""
  mtext(side=3, line=-1.1, at=2:6,text=sig, cex=0.8)
  mtext(side=3, text=paste("rho =",round(lsd$spearman[sp,"rho"],2),p2star(lsd$spearman[sp,"p.val"])), adj=1, cex= 0.7)
  title(main=sp, cex.main=0.8)
}

plot.lsd=function(lsd=LSDclassSR, var="SR", db=databp, sel.criteria = c( "spearman"),
                 spear = effect.summary) {
  
  par(mfrow=c(4,5), mar=c(2,2,2,2))
  a=0
  for (i in rownames(lsd$diff)){
    # only species with significant kruskal walis test :
    if (sel.criteria == "kruskal") sel <- lsd$kruskal.test[i,"p.val"] <=0.05
    if (sel.criteria == "anova") sel <- lsd$anova[i,"p.val"] <=0.05
    if (sel.criteria == "spear.Jecol") sel <- spear[i,"S.pval"] <=0.05
    if (sel.criteria == "spearman") sel <- lsd$spear[i,"p.val"] <=0.05 &  lsd$spear[i,"rho"] <0
    if (sel.criteria == "none") sel <- TRUE
    
    if (sel==T){
      plot.sp.lsd(sp= i, lsd=lsd, var=var, db=db)
      a=a+1
    }
  }
}



### Summarize results :
summary.lsd=function(M=LSDclassSR, data=species, group="ALIEN", graph=T) {
  data=data[rownames(M$diff),]
  if (is.null(group)) {
    sub=M
    G=1
  }
  
  if (!is.null(group)) { 
    sub=list()
    i=1
    G=length(unique(data[,group]))
    for (g in sort(unique(data[,group]))) {
      sub[[i]]=lapply(M, FUN= function(x) x[rownames(x) %in% rownames(data)[data[,group]==g],] )
      i=i+1
    }
    names=paste(group, ":",sort(unique(data[,group])), sep="")
  }                                 
  
  par(mfrow=c(G,2), oma=c(4,6,6,2), mar=c(2,2,1,1))
  
  out=NULL   
  for (j in 1:G) {
    S=sub[[j]]
    
    # number of potentieally signif species occurrence per class
    n.sp=sapply(names(S$diff), FUN= function(x) {
      length(which(!is.na(S$diff[,x])))
    })
    
    # number of significant negative effects per class
    n.impact=sapply(names(S$diff), FUN= function(x) {
      length(which(S$diff[,x] + S$leastdiff[,x]<0))
    })
    
    # number of significant negative effects per class
    n.positives=sapply(names(S$diff), FUN= function(x) {
      length(which(S$diff[,x] + S$leastdiff[,x]>0))
    })
    
    # proportion of signifif negative effects per class
    p.impact=n.impact / n.sp
    
    # number of times each class is the threshold
#     freq.thr= (f <- function() {
#       tmp=table(S$threshold$threshold)
#       y=as.vector(rep(0,5)) 
#       names(y)=names(S$diff)
#       y[as.numeric(names(tmp))-1]=tmp
#       return(y)
#     })()
    
    
    freq.thr =table(as.factor(S$threshold$threshold))[match(2:6,names(table(as.factor(S$threshold$threshold))))]
    names(freq.thr)=names(n.sp)
    freq.thr[is.na(freq.thr)]=0
    prop.thr = freq.thr / n.sp
    if(sum(freq.thr, na.rm=T)!=0) perc.thr = freq.thr / sum(freq.thr, na.rm=T)
    if(sum(freq.thr, na.rm=T)==0) perc.thr = freq.thr
    
    # proportion of times a class is a threshold
    prop.thr = freq.thr / colSums(!is.na(S$diff))
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

#test signif negative difference :
# tmp=LSDclass.test()
# (tmp$least + tmp$diff)<0
# table(rowSums((tmp$least + tmp$diff)<0, na.rm=T))  # number of species 
# # 0   1   2   3   4 
# # 201   3   4   6   1   in total only 14 species which have significant anova test 
# 
# LSDclassSRsig=names(which((rowSums((tmp$least + tmp$diff)<0, na.rm=T))>0))
# 
# # vs 25 species with a kruskal wallis + ttests
#  table(LSDclassSR$threshold$threshold)
# # 2 3 4 5 6 
# # 7 5 7 5 1 