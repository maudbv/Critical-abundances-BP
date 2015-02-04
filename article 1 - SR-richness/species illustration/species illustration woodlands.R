# REpresenting class-SR patterns for individual species in woodLAND sites

db=databp[databp$PlotName%in% woodlands,]

############ Negative IMPACTS ############
# Summary table of effects to be exported 
#write.csv(wood.effect.summary[wood.effect.summary$total.effect>0,], file='alleffects.csv')


# species with significant negative impact
negatives=na.omit(as.character(wood.effect.summary$SpeciesCode[wood.effect.summary$effect==-1]))

# Graph
x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in% negatives,]
tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(1,2), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  if(length(unique(x$domlevels))<2) plot(SR~as.factor(domlevels), x, border="grey")
  if(length(unique(x$domlevels))>2){
  co=cor.test(x$SR,x$domlevels, method='spearman')
 
  k=kruskal.test(x$SR,x$domlevels)
   
#   if (k$p.val<0.05 & co$est<0) plot(SR~as.factor(domlevels), x, col="forestgreen")
#   if (k$p.val<0.05 & co$est>0) plot(SR~as.factor(domlevels), x, col="firebrick")
#   if (k$p.val>0.05) plot(SR~as.factor(domlevels), x, col="grey")
  
  if (i %in% natives) plot(SR~as.factor(domlevels), x, col="forestgreen", varwidth=T)
  if (i %in% aliens) plot(SR~as.factor(domlevels), x, col="firebrick", varwidth=T)
  title(main=as.character(i))
  mtext(3,line=-1, adj=0.1,cex=0.6,
    text=paste('S:', round(wood.effect.summary$rho[wood.effect.summary$SpeciesCode==i],2),
    p2star(wood.effect.summary$S.pval[wood.effect.summary$SpeciesCode==i])))
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
    text=paste('K:',p2star(wood.effect.summary$k.pval[wood.effect.summary$SpeciesCode==i])))

  }
  }
  
# Graph
  x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in% negatives,]
tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(1,2), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
    x=tmp[as.character(tmp$SpeciesCode)==i,]
    plot(SR~domlevels, x, type='n', xlim=c(1,7))
    
    fit=loess(SR~domlevels, x, span=2)
    newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.1))
#     xx=c(t(newdata),newdata[dim(newdata)[1]:1,])
#     yy=c(predict(fit,newdata, se=T)$fit + 2*predict(fit,newdata, se=T)$se.fit ,
#     predict(fit,newdata, se=T)$fit[dim(newdata)[1]:1] -  2*predict(fit,newdata, se=T)$se.fit[dim(newdata)[1]:1] )
#     polygon(xx,yy,col='grey', border=NA)
    points(x$domlevels,x$SR,pch=20,  col=c("forestgreen", 'firebrick')[as.numeric(i %in% aliens)+1])
    
    lines(t(newdata),predict(fit,newdata, se=T)$fit)
    
    title(main=as.character(i))
    mtext(3,line=-1, adj=0.1,cex=0.6,
    text=paste('S:', round(wood.effect.summary$rho[wood.effect.summary$SpeciesCode==i],2),
    p2star(wood.effect.summary$S.pval[wood.effect.summary$SpeciesCode==i])))
    
    mtext(3,line=-1.7, adj=0.1,cex=0.6,
    text=paste('K:',p2star(wood.effect.summary$k.pval[wood.effect.summary$SpeciesCode==i])))
    
    }

## Impact on natives
  x11()
tmp=db
#tmp=tmp[tmp$SpeciesCode %in% negatives,]
tmp=tmp[tmp$SpeciesCode %in% as.character(wood.effect.summary$SpeciesCode[wood.effect.summary$nat.effect==-1 & wood.effect.summary$effect>-1]),]

tmp=tmp[order(-tmp$ALIEN),]
 par(mfrow=c(2,4), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRnat~domlevels, x, type='n', xlim=c(1,7))

  if (wood.effect.summary$nat.effect[wood.effect.summary$SpeciesCode==i]==-1) {
      fit=loess(SRnat~domlevels, x, span=2)
      newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.1))
#       xx=c(t(newdata),newdata[dim(newdata)[1]:1,])
#       yy=c(predict(fit,newdata, se=T)$fit + 2*predict(fit,newdata, se=T)$se.fit ,
#       predict(fit,newdata, se=T)$fit[dim(newdata)[1]:1] -  2*predict(fit,newdata, se=T)$se.fit[dim(newdata)[1]:1] ) 
#       polygon(xx,yy,col='grey', border=NA)
      lines(t(newdata),predict(fit,newdata, se=T)$fit)
      }
 points(x$domlevels,x$SRnat,pch=20,  col=c("forestgreen", 'firebrick')[as.numeric(i %in% aliens)+1])

  title(main=as.character(i))
  mtext(3,line=-1, adj=0.1,cex=0.6,
        text=paste('S:', round(wood.effect.summary$nat.rho[wood.effect.summary$SpeciesCode==i],2),
                   p2star(wood.effect.summary$nat.S.pval[wood.effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(wood.effect.summary$nat.k.pval[wood.effect.summary$SpeciesCode==i])))
                 
  }
 
## Impact on aliens
  x11()
tmp=db
#tmp=tmp[tmp$SpeciesCode %in% negatives,]
tmp=tmp[tmp$SpeciesCode %in% as.character(wood.effect.summary$SpeciesCode[wood.effect.summary$ali.effect==-1 & wood.effect.summary$effect>-1]),]

tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(1,4), mar=c(3,2,2,1))
#par(mfrow=c(1,2), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRali~domlevels, x, type='n', xlim=c(1,7))

  if (wood.effect.summary$ali.effect[wood.effect.summary$SpeciesCode==i]==-1) {
      fit=loess(SRali~domlevels, x, span=2)
      newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.1))
#       xx=c(t(newdata),newdata[dim(newdata)[1]:1,])
#       yy=c(predict(fit,newdata, se=T)$fit + 2*predict(fit,newdata, se=T)$se.fit ,
#       predict(fit,newdata, se=T)$fit[dim(newdata)[1]:1] -  2*predict(fit,newdata, se=T)$se.fit[dim(newdata)[1]:1] ) 
#       polygon(xx,yy,col='grey', border=NA)
      lines(t(newdata),predict(fit,newdata, se=T)$fit)
      }
 points(x$domlevels,x$SRali,pch=20,  col=c("forestgreen", 'firebrick')[as.numeric(i %in% aliens)+1])

  title(main=as.character(i))
  mtext(3,line=-1, adj=0.1,cex=0.6,
        text=paste('S:', round(wood.effect.summary$ali.rho[wood.effect.summary$SpeciesCode==i],2),
                   p2star(wood.effect.summary$ali.S.pval[wood.effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(wood.effect.summary$ali.k.pval[wood.effect.summary$SpeciesCode==i])))
                 
  }
 
############ POSITIVE IMPACTS ############

# species with significant positive impact


  ## Impact on natives
  x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in% as.character(wood.effect.summary$SpeciesCode[wood.effect.summary$nat.effect==1]),]
tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(2,3), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRnat~domlevels, x, type='n', xlim=c(1,7))

  if (wood.effect.summary$nat.effect[wood.effect.summary$SpeciesCode==i]==1) {
      fit=loess(SRnat~domlevels, x, span=2)
      newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.1))
#       xx=c(t(newdata),newdata[dim(newdata)[1]:1,])
#       yy=c(predict(fit,newdata, se=T)$fit + 2*predict(fit,newdata, se=T)$se.fit ,
#       predict(fit,newdata, se=T)$fit[dim(newdata)[1]:1] -  2*predict(fit,newdata, se=T)$se.fit[dim(newdata)[1]:1] ) 
#       polygon(xx,yy,col='grey', border=NA)
      lines(t(newdata),predict(fit,newdata, se=T)$fit)
      }
 points(x$domlevels,x$SRnat,pch=20,  col=c("forestgreen", 'firebrick')[as.numeric(i %in% aliens)+1])

  title(main=as.character(i))
  mtext(3,line=-1, adj=0.1,cex=0.6,
        text=paste('S:', round(wood.effect.summary$nat.rho[wood.effect.summary$SpeciesCode==i],2),
                   p2star(wood.effect.summary$nat.S.pval[wood.effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(wood.effect.summary$nat.k.pval[wood.effect.summary$SpeciesCode==i])))
                 
  }
 
## Impact on aliens
  x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in% as.character(wood.effect.summary$SpeciesCode[wood.effect.summary$ali.effect==1]),]

tmp=tmp[order(-tmp$ALIEN),]
 par(mfrow=c(3,5), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRali~domlevels, x, type='n', xlim=c(1,7))

  if (wood.effect.summary$ali.effect[wood.effect.summary$SpeciesCode==i]==1) {
      fit=loess(SRali~domlevels, x, span=2)
      newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.1))
#       xx=c(t(newdata),newdata[dim(newdata)[1]:1,])
#       yy=c(predict(fit,newdata, se=T)$fit + 2*predict(fit,newdata, se=T)$se.fit ,
#       predict(fit,newdata, se=T)$fit[dim(newdata)[1]:1] -  2*predict(fit,newdata, se=T)$se.fit[dim(newdata)[1]:1] ) 
#       polygon(xx,yy,col='grey', border=NA)
      lines(t(newdata),predict(fit,newdata, se=T)$fit)
      }
 points(x$domlevels,x$SRali,pch=20,  col=c("forestgreen", 'firebrick')[as.numeric(i %in% aliens)+1])

  title(main=as.character(i))
  mtext(3,line=-1, adj=0.1,cex=0.6,
        text=paste('S:', round(wood.effect.summary$ali.rho[wood.effect.summary$SpeciesCode==i],2),
                   p2star(wood.effect.summary$ali.S.pval[wood.effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(wood.effect.summary$ali.k.pval[wood.effect.summary$SpeciesCode==i])))
                 
  }