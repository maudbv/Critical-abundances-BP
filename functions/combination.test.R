
simul1.combinations=function(db=databp,status='ALIEN', out=c('rhos', 'pval'),envplot=envplot,zeros=F, min=10,alpha=0.01, nreps=99)  {

require(doBy)

  db$Sp.occurence=table(as.character(db$SpeciesCode))[match(as.character(db$SpeciesCode),
                                                            unique(as.character(db$SpeciesCode)))]
  db=db[db$Sp.occurence>min,]

nat.rhos=data.frame(matrix(NA, ncol=nreps+1, nrow=length(unique(db$SpeciesCode))))
ali.rhos=data.frame(matrix(NA, ncol=nreps+1, nrow=length(unique(db$SpeciesCode))))

# OBSERVED values

    classSR.nat=classSR.fun(db=db, var="SRnat",envplot=envplot,zeros=zeros, min.occur=min, min.class=3, alpha=alpha)
    classSR.ali=classSR.fun(db=db, var="SRali", envplot=envplot,zeros=zeros, min.occur=min, min.class=3, alpha=alpha)

    nat.rhos[,1]= classSR.nat$rho
    ali.rhos[,1]= classSR.ali$rho

## SIMULATIONS  
  for (i in 1:nreps) {
    
    #Null shuffling of abundance classes within sites
    nul=tapply(as.character(db$domlevels),as.character(db$PlotName), FUN=function(x) sample(x))
    nul.db=db
    nul.db$domlevels=as.numeric(as.character(unlist(nul)))
  
    tmp=nul.db
    null.classSR.nat=classSR.fun(db=tmp, var="SRnat", envplot=envplot,zeros=zeros, min.occur=min, min.class=3, alpha=alpha)
    nat.rhos[,i+1]=  null.classSR.nat$rho
    null.classSR.ali=classSR.fun(db=tmp, var="SRali", envplot=envplot,zeros=zeros, min.occur=min, min.class=3, alpha=alpha)
    ali.rhos[,i+1]=  null.classSR.ali$rho
}   

  rownames(nat.rhos)=   rownames(ali.rhos)= unique(db$SpeciesCode) 
  
  return(list(nat.rhos,ali.rhos))
}

############ Plotting results

results.sim1.combination=function(sim=sim1.combi, by.status=T, filter=T) {
  x=-sim[[1]]
  y=-sim[[2]]
  
if(filter==T){  
  tmp=which(abs(x[,1])<0.05 | abs(y[,1])<0.05)
  x[tmp,1]=NA
  y[tmp,1]=NA
}
  
  if (by.status==F) {
  # compare frequency of combinations
obs.freq=as.data.frame(table(neg.NR=x[,1]<=0, neg.AR=y[,1]<=0))
null.freq=NULL
for (i in 2:dim(x)[2]){
  tmp=as.data.frame(table(neg.NR=x[,i]<=0, neg.AR=y[,i]<=0))
null.freq=cbind(null.freq,tmp$Freq)
}
es=(rowSums((null.freq[,] < obs.freq[,4])) + rowSums((null.freq[,] == obs.freq[,3]))/2)/length(null.freq[1,])
p=(1-abs(2*(es-0.5)))

  # graph
  plot(c(0.5,4.5),c(min(obs.freq[,3])-5,max(obs.freq[,3])+10), type='n', ann=F, axes=F)
  quant=apply(null.freq,1, quantile, probs=c(0.0275,0.9725))
  arrows(1:4, quant[1,],1:4, quant[2,], col='grey50',angle = 90,code =3, length=0.1)  
  boxplot(t(null.freq), pch='.',pars = list(boxwex = 0.5, staplewex = 0, outwex = 0),add=T,
            lty='solid',border='grey50',range=0.1, outline=F, 
            names=c('all positive', '- NR\n+ AR', '+ NR\n- AR', 'all negative'), cex.axis=0.7)
  points(1:4,obs.freq[,3], pch=20, cex=2)
  mtext(side=3, text=paste('P=',round(p,4)), at=1:4, cex=0.7)
  # mtext(side=1,line=3, font=3, text=c("Facilitation\nOptimal habitat", 'Native-specific\nnegative impact?', 'Resistance','Competition'), 
  #       at=1:4, cex=0.7, col='grey40')
}

if (by.status==T) {
  #BY STATUS compare frequency of combinations
obs.freq=as.data.frame(table(neg.NR=x[,1]<=0, neg.AR=y[,1]<=0,alien=rownames(x)%in%aliens))
null.freq=NULL
for (i in 2:dim(x)[2]){
  tmp=as.data.frame(table(neg.NR=x[,i]<=0, neg.AR=y[,i]<=0,alien=rownames(x)%in%aliens))
null.freq=cbind(null.freq,tmp$Freq)
}
es=(rowSums((null.freq[,] < obs.freq[,4])) + rowSums((null.freq[,] == obs.freq[,3]))/2)/length(null.freq[1,])
p=(1-abs(2*(es-0.5)))

  #GRaph
  plot(c(0.5,8.5),c(min(obs.freq[,4])-5,max(obs.freq[,4])+10), type='n', ann=F, axes=F, bty='o')
  quant=apply(null.freq,1, quantile, probs=c(0.0275,0.9725))
  arrows(1:8, quant[1,],1:8, quant[2,], col='grey50',angle = 90,code =3, length=0.1)  
  boxplot(t(null.freq), pch='.',pars = list(boxwex = 0.5, staplewex = 0, outwex = 0),add=T,
            lty='solid',border='grey50',range=0.1, outline=F, col=c(rep('white',4),rep('grey',4)),
            names=NULL,xaxt='n', cex.axis=0.8, las=1)
  axis(side=1,at=1:8,label=rep(c('all\npositive', '- NR\n+ AR', '+ NR\n- AR', 'all\nnegative'),2),
       tcl=-0.3,cex=0.7, padj=0.5)
  points(1:8,obs.freq[,4], pch=20, cex=2)
  # mtext(side=3, text=paste('P=',round(p,4)), at=1:8, cex=0.7)
  mtext(side=3, text=sapply(p, p2star), at=1:8, cex=0.7, line=-1)
  abline(v=4.5, col="grey", lty="dotted")  
  mtext(side=3, text=c("Natives",  "Aliens"), at=c(2.5,6.5), cex=1)
}
  
  
# Table of results
output=cbind(obs.freq, mean.null=apply(null.freq, 1, mean, na.rm=T),sd.null=apply(null.freq, 1,sd, na.rm=T),t(quant), ES=es, p)
return(output)
  
  
#   # compare correlations
#   
#   x0=apply(x, 1, function(t) {
#     p=NA
#     if (!is.na(t[1])) {
#     p=(sum(t< as.numeric(t[1]), na.rm=T) + sum(t==t[1], na.rm=T)/2)/length(na.omit(t))
#     }
#     return(p)
#   })
#   
#   
#   y0=apply(y, 1, function(t) {
#     p=NA
#     if (!is.na(t[1])) {
#     p=(sum(t< as.numeric(t[1]), na.rm=T) + sum(t==t[1], na.rm=T)/2)/length(na.omit(t))
#     }
#     return(p)
#   })
#  
#   plot( y0~ x0)
#   obs.fit=summary(lm(y[,1]~ x[,1]))
#   combi.corr=rbind(r=obs.fit$adj.r.squared, p=obs.fit$coef[2,4])
# null.freq=NULL
# for (i in 2:dim(x)[2]){
#  null.fit=summary(lm( y0[,i]~ x0[,i]))
#  combi.corr=cbind(combi.corr, c(null.fit$adj.r.squared,null.fit$coef[2,4]))
# }
#   
#   hist(combi.corr[1,])
#   abline(v=combi.corr[1,1], col='red')
#   abline(v=quantile(combi.corr[1,], 0.95))
#  
#   hist(combi.corr[2,])
#   abline(v=combi.corr[2,2], col='red')
#   abline(v=quantile(combi.corr[2,], 0.05))

  
}


