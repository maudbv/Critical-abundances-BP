# REpresenting class-SR patterns for individual species in GRASSLAND sites

db=databp[databp$PlotName%in% grasslands,]
db$Sp.occurence.new=table(as.character(db$SpeciesCode))[match(db$SpeciesCode, unique(db$SpeciesCode))]
############ PODHAL
par(mfrow=c(2,3), oma=c(3,3,1,1), mar=c(3,2,2,1), las=1)
tmp=db
tmp$domlevels= - tmp$domlevels
i='PODHAL'
x=tmp[as.character(tmp$SpeciesCode)==i,]
    t=rep(1, length(x$SR))
    t[x$vegtype=='W'& !is.na(x$vegtype)]=2
    t[x$vegtype=='G'& !is.na(x$vegtype)]=3
    p=21
    c=c('grey','forestgreen','goldenrod')[t]

    plot(SR~as.numeric(domlevels), x, type='n', xlim=c(-6,-1), ylim=c(0,65), xaxt='n')
    axis(1, at=-6:-1 , label=1:6)
    fit=lm(SR~domlevels, x)
    newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.01))
    points(x$domlevels,jitter(x$SR),col=c, pch=p)
    lines(t(newdata),predict(fit,newdata, se=T)$fit)
    title(main="Total richness")
    mtext(3,line=-1, adj=0.1,cex=0.6,
    text=paste('S:', -round(effect.summary$rho[grass.effect.summary$SpeciesCode==i],2),
    p2star(effect.summary$S.pval[grass.effect.summary$SpeciesCode==i])))

    plot(SRnat~as.numeric(domlevels), x, type='n', xlim=c(-6,-1), ylim=c(0,65), xaxt='n')
    axis(1, at=-6:-1 , label=1:6)
 fit=lm(SRnat~domlevels, x)
    newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.01))
    points(x$domlevels,jitter(x$SRnat),col=c, pch=p)
    lines(t(newdata),predict(fit,newdata, se=T)$fit) 
    title(main="Native richness")
    mtext(3,line=-1, adj=0.1,cex=0.6,
    text=paste('S:', -round(effect.summary$nat.rho[effect.summary$SpeciesCode==i],2),
    p2star(effect.summary$nat.S.pval[effect.summary$SpeciesCode==i])))

  plot(SRali~as.numeric(domlevels), x, type='n', xlim=c(-6,-1), ylim=c(0,65), xaxt='n')
    axis(1, at=-6:-1 , label=1:6)
 fit=lm(SRali~domlevels, x, span=2)
    points(x$domlevels,jitter(x$SRali),col=c, pch=p)
    title(main="Alien richness")
    mtext(3,line=-1, adj=0.1,cex=0.6,
    text=paste('S:', -round(effect.summary$ali.rho[effect.summary$SpeciesCode==i],2),
    p2star(effect.summary$ali.S.pval[effect.summary$SpeciesCode==i])))
 
############ Negative IMPACTS ############
# Summary table of effects to be exported 
#write.csv(grass.effect.summary[grass.effect.summary$total.effect>0,], file='alleffects.csv')


# species with significant negative impact
negatives=na.omit(as.character(grass.effect.summary$SpeciesCode[grass.effect.summary$effect==-1]))

# Graph
x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in% negatives,]
tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(2,5), mar=c(3,2,2,1))
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
    text=paste('S:', round(grass.effect.summary$rho[grass.effect.summary$SpeciesCode==i],2),
    p2star(grass.effect.summary$S.pval[grass.effect.summary$SpeciesCode==i])))
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
    text=paste('K:',p2star(grass.effect.summary$k.pval[grass.effect.summary$SpeciesCode==i])))

  }
  }
  
# Graph
  x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in% negatives,]
tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(2,5), mar=c(3,2,2,1))
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
    text=paste('S:', round(grass.effect.summary$rho[grass.effect.summary$SpeciesCode==i],2),
    p2star(grass.effect.summary$S.pval[grass.effect.summary$SpeciesCode==i])))
    
    mtext(3,line=-1.7, adj=0.1,cex=0.6,
    text=paste('K:',p2star(grass.effect.summary$k.pval[grass.effect.summary$SpeciesCode==i])))
    
    }

## Impact on natives
  x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in% negatives,]
tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(2,5), mar=c(3,2,2,1))
# par(mfrow=c(1,3), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRnat~domlevels, x, type='n', xlim=c(1,7))

  if (grass.effect.summary$nat.effect[grass.effect.summary$SpeciesCode==i]==-1) {
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
        text=paste('S:', round(grass.effect.summary$nat.rho[grass.effect.summary$SpeciesCode==i],2),
                   p2star(grass.effect.summary$nat.S.pval[grass.effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(grass.effect.summary$nat.k.pval[grass.effect.summary$SpeciesCode==i])))
                 
  }
 
## Impact on aliens
  x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in% negatives,]
tmp=tmp[tmp$SpeciesCode %in% as.character(grass.effect.summary$SpeciesCode[grass.effect.summary$ali.effect==-1 & grass.effect.summary$effect>-1]),]

tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(2,5), mar=c(3,2,2,1))
# par(mfrow=c(1,3), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRali~domlevels, x, type='n', xlim=c(1,7))

  if (grass.effect.summary$ali.effect[grass.effect.summary$SpeciesCode==i]==-1) {
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
        text=paste('S:', round(grass.effect.summary$ali.rho[grass.effect.summary$SpeciesCode==i],2),
                   p2star(grass.effect.summary$ali.S.pval[grass.effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(grass.effect.summary$ali.k.pval[grass.effect.summary$SpeciesCode==i])))
                 
  }
 
############ POSITIVE IMPACTS ############

# species with significant positive impact
positives=as.character(grass.effect.summary$SpeciesCode[grass.effect.summary$effect==1])

# Graph
tmp=db
tmp=tmp[tmp$SpeciesCode %in% positives,]
par(mfrow=c(4,5), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  if(length(unique(x$domlevels))<2) plot(SR~as.factor(domlevels), x, border="grey")
  if(length(unique(x$domlevels))>2){
  co=cor.test(x$SR,x$domlevels, method='spearman')
 
  k=kruskal.test(x$SR,x$domlevels)
   
#   if (k$p.val<0.05 & co$est<0) plot(SR~as.factor(domlevels), x, col="forestgreen")
#   if (k$p.val<0.05 & co$est>0) plot(SR~as.factor(domlevels), x, col="firebrick")
#   if (k$p.val>0.05) plot(SR~as.factor(domlevels), x, col="grey")
  
  if (i %in% natives) plot(SR~as.factor(domlevels), x, col="forestgreen")
  if (i %in% aliens) plot(SR~as.factor(domlevels), x, col="firebrick")
  title(main=as.character(i))
  mtext(3,text=round(co$est,2), adj=1,cex=0.6)
  }
  }



# Graph
  x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in%positives,]
tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(1,5), mar=c(3,2,2,1))
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
    text=paste('S:', round(grass.effect.summary$rho[grass.effect.summary$SpeciesCode==i],2),
    p2star(grass.effect.summary$S.pval[grass.effect.summary$SpeciesCode==i])))
    
    mtext(3,line=-1.7, adj=0.1,cex=0.6,
    text=paste('K:',p2star(grass.effect.summary$k.pval[grass.effect.summary$SpeciesCode==i])))
    
    }

  ## Impact on natives
  x11()
tmp=db
tmp=tmp[tmp$SpeciesCode %in% positives,]

tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(1,3), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRnat~domlevels, x, type='n', xlim=c(1,7))

  if (grass.effect.summary$nat.effect[grass.effect.summary$SpeciesCode==i]==1) {
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
        text=paste('S:', round(grass.effect.summary$nat.rho[grass.effect.summary$SpeciesCode==i],2),
                   p2star(grass.effect.summary$nat.S.pval[grass.effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(grass.effect.summary$nat.k.pval[grass.effect.summary$SpeciesCode==i])))
                 
  }
 
## Impact on aliens
  x11()
tmp=db
tmp=tmp[na.omit(tmp$SpeciesCode %in% positives),]
#tmp=tmp[tmp$SpeciesCode %in% as.character(grass.effect.summary$SpeciesCode[grass.effect.summary$ali.effect==1 & grass.effect.summary$effect<1]),]

tmp=tmp[order(-tmp$ALIEN),]
 par(mfrow=c(1,5), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRali~domlevels, x, type='n', xlim=c(1,7))

  if (grass.effect.summary$ali.effect[grass.effect.summary$SpeciesCode==i]==1) {
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
        text=paste('S:', round(grass.effect.summary$ali.rho[grass.effect.summary$SpeciesCode==i],2),
                   p2star(grass.effect.summary$ali.S.pval[grass.effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(grass.effect.summary$ali.k.pval[grass.effect.summary$SpeciesCode==i])))
                 
  }