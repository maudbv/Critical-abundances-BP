classSR.fun=function(db=databp, var='SR',domclass='domlevels',zeros=F,env=envplot, min.occur=10,  min.class=3, alpha=0.01) {
  #  function testing all relationships between dominance leveles (class) and a variable of the community (e.g. SR)
  classSR=NULL
  db$Sp.occurence=table(as.character(db$SpeciesCode))[match(as.character(db$SpeciesCode),unique(as.character(db$SpeciesCode)))]
  tmp=db[db$Sp.occurence>min.occur,]
  
  for (i in unique(tmp$SpeciesCode) ) {
    x=tmp[as.character(tmp$SpeciesCode)==i,c(domclass,var,'PlotName')]
    names(x)=c('abun','var','PlotName')
    z= cbind(abun=NA,env[which(!env$PLOTID %in% x$PlotName),c(var,"PLOTID")])
    x=apply(x,2,as.numeric)
    z=apply(z,2,as.numeric)
    x=as.data.frame(rbind(x,z))

    x$abun=7-x$abun
    if (zeros) x$abun[is.na(x$abun)]=0
    
    nabun=length(na.omit(unique(x$abun[x$abun>0])))
  
    # if the species is not represented in at least 'min.abun' different levels of abundance :
    if(length(unique(x$abun)) < min.class) {
      r=c(unique(tmp$Sp.occurence[i]),nabun,
      max(x$abun, na.rm=T), min(x$abun[x$abun>0], na.rm=T),
      NA,NA, NA,NA,NA)
      }
#       print(length(unique(x$abun)))
    
    # if the species is represented in sufficient different levels of abundance :
    if(length(unique(x$abun))>=min.class) {
      co=cor.test(x$var,x$abun, method='spearman', exact=F)
      k=kruskal.test(x$var,x$abun)  
      r=c(unique(tmp$Sp.occurence[i]),nabun,
      max(x$abun, na.rm=T), min(x$abun[x$abun>0], na.rm=T),
      k$stat, k$param, k$p.value,co$estim, co$p.val)
      }
    classSR=rbind(classSR,r)
    }
  classSR=data.frame(as.character(unique(tmp$SpeciesCode)),classSR, row.names=as.character(unique(tmp$SpeciesCode)))
  names(classSR)=c('SpeciesCode','Sp.occurence','nb.class','max','min','K.chi2','df','k.pval','rho','S.pval')
  
  classSR=merge(classSR, species,by.x='SpeciesCode', by.y='Sp.code', all.x=T)
  classSR$effect=0
  #classSR$effect[which(classSR$k.pval<=alpha | classSR$S.pval<=alpha)]=1
  classSR$effect[which(classSR$S.pval<=alpha)]=1
  classSR$effect[which(classSR$effect==1 & classSR$rho<0)]=-1
  classSR$effect[is.na(classSR$rho)]=NA
  return(classSR)
}

