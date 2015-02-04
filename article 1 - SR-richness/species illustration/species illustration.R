# REpresenting class-SR patterns for individual species


## frequency of different abundance classes
classfreq=rbind(table(as.numeric(databp$domlevels)),
                table(as.numeric(databp$domlevels[databp$ALIEN==1])),
                table(as.numeric(databp$domlevels[databp$ALIEN==0])))
classfreq=classfreq/rowSums(classfreq)
barplot(classfreq[2:3,], beside=T, )
summary(table(databp$domlevels,databp$ALIEN))

# 
# classfreq=as.data.frame(table(databp$SpeciesCode,databp$domlevels))
# classfreq$ALIEN=species$ALIEN[match(rownames(classfreq), species$Sp.code)]
# classfreq[,1:6]=classfreq[,1:6]/rowSums(classfreq[,1:6])
# mean.freq=data.frame(natives=colMeans(classfreq[classfreq$ALIEN==0,1:6]),
#                      aliens=colMeans(classfreq[classfreq$ALIEN==1,1:6], na.rm=T))
# barx=barplot(t(as.matrix(mean.freq)), beside=T)
# sd.freq=data.frame(natives=apply(classfreq[classfreq$ALIEN==0,1:6], 2, sd, na.rm=T),
#                      aliens=apply(classfreq[classfreq$ALIEN==1,1:6], 2, sd, na.rm=T))
# segments(barx,t(as.matrix(mean.freq)) + t(as.matrix(sd.freq)),barx ,t(as.matrix(mean.freq)))


############ plotting FUNCTION
plot.illustration=function(i='LOLPER',var="SR", db=envplot, X=comm, zeros=F, 
                           c="grey40",xaxis=T, yaxis=T, vegtype=F) {
  x=data.frame(A=X[match(as.character(db$PLOTID),rownames(X)),i],var=na.omit(as.vector(db[,var])))
  x$A=7-x$A
  if (zeros)  x$A[x$A==7]=0
  if (!zeros)  x$A[x$A==7]=NA
      
  ## distinguish vegetation type in colors
  p=21
    if (vegtype) {
        t=rep(1, length(x[,var]))
        t[x$vegtype=='W'& !is.na(x$vegtype)]=2
        t[x$vegtype=='G'& !is.na(x$vegtype)]=3
        c=c('grey','black','grey50')[t]
        bg=c(NA,'black','grey')[t]
      }
    else {
      bg=NA
      c=c
    }
   
 ## plotting
    plot(var ~ A, data=x, type='n', xlim=c(0,6), ylim=c(0,65), xaxt='n', yaxt='n', bty='l')
#     if (xaxis) axis(1, at=1:6 , label=1:6, tcl=0.3, mgp=c(3,0.1,0),cex.axis=0.9)
#     if (yaxis) axis(2, las=1, cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
    if (xaxis==F) axis(1, at=1:6 , label=NA, tcl=0.3, mgp=c(3,0.1,0),cex.axis=0.9)
    if (yaxis==F) axis(2, las=1, label=NA,cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
  
    if (xaxis==T) axis(1, at=1:6 , label=1:6, tcl=0.3, mgp=c(3,0.1,0),cex.axis=0.9)
    if (yaxis==T) axis(2, las=1, cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
  
    fit=lm(var ~ A, data=x)
    newdata=data.frame(domlevels=seq(min(x$A, na.rm=T),max(x$A, na.rm=T),0.01))
    points(x$A,jitter(x$var),col=c, pch=p,bg=bg)
 
## spearman rank test correlation 
    co=cor.test(x$A, x$var, method="spearman")
    mtext(3,line=-1.1, adj=0.97,cex=0.7,
    text=substitute(italic(rho)*' = '*x, list(x=paste(round(co$estim,2),p2star(co$p.val), sep=''))))
   if(co$p.val<0.01) lines(na.omit(x$A),predict(fit))
  }

# x11()
# par(mfcol=c(3,4), oma=c(3,10,4,1), mar=c(1,1,0,0), las=1)
# plot.illustration(i='LOLPER',var="SR", xaxis=F, yaxis=T)
# plot.illustration(i='LOLPER',var="SRnat", xaxis=F, yaxis=T)
# plot.illustration(i='LOLPER',var="SRali", xaxis=T, yaxis=T)
#  
# plot.illustration(i='DIGPUR',var="SR", xaxis=F, yaxis=F)
# plot.illustration(i='DIGPUR',var="SRnat", xaxis=F, yaxis=F)
# plot.illustration(i='DIGPUR',var="SRali", xaxis=T, yaxis=F)
# 
# plot.illustration(i='MELRAM',var="SR", xaxis=F, yaxis=F)
# plot.illustration(i='MELRAM',var="SRnat", xaxis=F, yaxis=F)
# plot.illustration(i='MELRAM',var="SRali", xaxis=T, yaxis=F)
# 
# plot.illustration(i='COPRHA',var="SR", xaxis=F, yaxis=F)
# plot.illustration(i='COPRHA',var="SRnat", xaxis=F, yaxis=F)
# plot.illustration(i='COPRHA',var="SRali", xaxis=T, yaxis=F)
# # 
# mtext(1,text="Abundance class", outer=T, cex=0.8)      
# mtext(2,text=expression(italic("L. perenne\n(Alien grass)")), outer=T, at=0.85, line=4, cex=0.8, las=1, adj=0.5)
# mtext(2,text=expression(italic("M. ramiflorus\n(Native tree)")), outer=T, at=0.35, line=4, cex=0.8, las=1, adj=0.5)
# mtext(2,text="Number of species", outer=T, at=0.52, line=9, cex=0.9, las=0)
# 
# mtext(3,text="Total richness", outer=T, at=0.18, line=0, cex=0.8, las=0)
# mtext(3,text="Native richness", outer=T, at=0.51, line=0, cex=0.8, las=0)
# mtext(3,text="Alien richness", outer=T, at=0.85, line=0, cex=0.8, las=0)
# 
# mtext(1,text="Abundance level", outer=T, line=0.5,cex=0.8)      
# mtext(3,text=expression(italic("L. perenne\n(Alien grass)")), outer=T, at=0.1, line=1, cex=0.8, las=1, adj=0.5)
# mtext(3,text=expression(italic("D. purpurea\n(Alien herb)")), outer=T, at=0.35, line=1, cex=0.8, las=1, adj=0.5)
# mtext(3,text=expression(italic("M. ramiflorus\n(Native tree)")), outer=T, at=0.6, line=1, cex=0.8, las=1, adj=0.5)
# mtext(3,text=expression(italic("C. rotundifolia\n(Native shrub)")), outer=T, at=0.85, line=1, cex=0.8, las=1, adj=0.5)
# 
# mtext(2,text="Species richness", outer=T, at=0.52, line=8, cex=0.9, las=0)
# 
# mtext(2,text="Alien\nrichness", outer=T, at=0.18, line=2, cex=0.8, las=1)
# mtext(2,text="Native\nrichness", outer=T, at=0.51, line=2, cex=0.8, las=1)
# mtext(2,text="Total\nrichness", outer=T, at=0.85, line=2, cex=0.8, las=1)

#New version  
x11()
par(mfcol=c(3,4), oma=c(3,10,4,1), mar=c(1,1,0,0), las=1)
plot.illustration(i='CRIMUR',var="SR", xaxis=F, yaxis=T, c="black")
plot.illustration(i='CRIMUR',var="SRnat", xaxis=F, yaxis=T, c="black")
plot.illustration(i='CRIMUR',var="SRali", xaxis=T, yaxis=T, c="black")
 
plot.illustration(i='RUMACE',var="SR", xaxis=F, yaxis=F, c="black")
plot.illustration(i='RUMACE',var="SRnat", xaxis=F, yaxis=F, c="black")
plot.illustration(i='RUMACE',var="SRali", xaxis=T, yaxis=F, c="black")

plot.illustration(i='GRILIT',var="SR", xaxis=F, yaxis=F, c="black")
plot.illustration(i='GRILIT',var="SRnat", xaxis=F, yaxis=F, c="black")
plot.illustration(i='GRILIT',var="SRali", xaxis=T, yaxis=F, c="black")

plot.illustration(i='COPROT',var="SR", xaxis=F, yaxis=F, c="black")
plot.illustration(i='COPROT',var="SRnat", xaxis=F, yaxis=F, c="black")
plot.illustration(i='COPROT',var="SRali", xaxis=T, yaxis=F, c="black")

mtext(1,text="Abundance class", outer=T, line=0.5,cex=0.8)      
mtext(3,text=expression(italic("a) C. murinum\n   (Alien grass)")), outer=T, at=0.12, line=1, cex=0.8, las=1, adj=0.5)
mtext(3,text=expression(italic("b) R. acetosella\n   (Alien herb)")), outer=T, at=0.37, line=1, cex=0.8, las=1, adj=0.5)
mtext(3,text=expression(italic("c) G. littoralis\n   (Native tree)")), outer=T, at=0.62, line=1, cex=0.8, las=1, adj=0.5)
mtext(3,text=expression(italic("d) C. rotundifolia\n   (Native shrub)")), outer=T, at=0.87, line=1, cex=0.8, las=1, adj=0.5)

mtext(2,text="Species richness", outer=T, at=0.52, line=8, cex=0.9, las=0)

mtext(2,text="Alien\nrichness", outer=T, at=0.18, line=2, cex=0.8, las=1)
mtext(2,text="Native\nrichness", outer=T, at=0.51, line=2, cex=0.8, las=1)
mtext(2,text="Total\nrichness", outer=T, at=0.85, line=2, cex=0.8, las=1)





# With zeros

plot.illustration =function(i='LOLPER',var="SR",db=envplot, X=comm) {
#    db=db[db$PLOTID %in% woodlands,]
  x=data.frame(A=X[match(as.character(db$PLOTID),rownames(X)),i],var=na.omit(db[,var]))
  x$A=7-x$A
  x$A[x$A==7]=0
    plot(var ~ A, data=x, type='n',ylim=c(0,65), xlim=c(0,7),xaxt='n', yaxt='n', bty='l')
    axis(1, at=0:6 , label=0:6, tcl=0.3, mgp=c(3,0.3,0),cex.axis=1,)
    axis(2, las=1, cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
    fit=lm(var ~ A, data=x)
    newdata=data.frame(domlevels=seq(min(x$A, na.rm=T),max(x$A, na.rm=T),0.01))
    points(x$A,jitter(x$var))
  
   co=cor.test(x$A, x$var, method="spearman")
    mtext(3,line=-1, adj=0.95,cex=0.6,
    text=substitute(italic(rho)*' = '*x, list(x=paste(round(co$estim,2),p2star(co$p.val), sep=''))))
   if(co$p.val<0.01) lines(na.omit(x$A),predict(fit))
}
  
# x11()
par(mfrow=c(2,3), oma=c(3,10,1,1), mar=c(2,2,1,1), las=1)
plot.illustration(i='LOLPER',var="SR")
plot.illustration(i='LOLPER',var="SRnat")
plot.illustration(i='LOLPER',var="SRali")
 
plot.illustration(i='MELRAM',var="SR")
plot.illustration(i='MELRAM',var="SRnat")
plot.illustration(i='MELRAM',var="SRali")

mtext(1,text="Abundance class", outer=T, cex=0.8)      
mtext(2,text=expression(italic("L. perenne\n(Alien grass)")), outer=T, at=0.75, line=5, cex=0.8, las=1, adj=0.5)
mtext(2,text=expression(italic("M. ramiflorus\n(Native tree)")), outer=T, at=0.25, line=5, cex=0.8, las=1, adj=0.5)
mtext(2,text="Number of species", outer=T, at=0.52, line=0, cex=0.8, las=0)

mtext(3,text="Total richness", outer=T, at=0.18, line=0, cex=0.8, las=0)
mtext(3,text="Native richness", outer=T, at=0.51, line=0, cex=0.8, las=0)
mtext(3,text="Alien richness", outer=T, at=0.85, line=0, cex=0.8, las=0)
          
 
############ Lolium perenne vs MELRAM : effects on SRwoody and SR legumes

par(mfrow=c(2,2), oma=c(3,3,1,1), mar=c(3,2,2,1), las=1)
tmp=databp
i='LOLPER'
x=tmp[as.character(tmp$SpeciesCode)==i,]
    t=rep(1, length(x$SR))
    t[x$vegtype=='W'& !is.na(x$vegtype)]=2
    t[x$vegtype=='G'& !is.na(x$vegtype)]=3
    p=21
    c=c('grey','forestgreen','goldenrod')[t]

    plot(SR.wood~as.numeric(domlevels), x, type='n', xlim=c(1,6))
    fit=lm(SR.wood~domlevels, x)
    newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.01))
    points(x$domlevels,jitter(x$SR.wood),col=c, pch=p)
    lines(t(newdata),predict(fit,newdata, se=T)$fit)
    title(main="Woody richness")
    mtext(3,line=-1, adj=0.1,cex=0.6,
    text=paste('S:', round(cor.test(x$domlevels,x$SR.wood, method='spearman')$estimate,2),
    p2star(cor.test(x$domlevels,x$SR.wood, method='spearman')$p.val)))

    plot(SR.legumes~as.numeric(domlevels), x, type='n', xlim=c(1,6))
    fit=lm(SR.legumes~domlevels, x)
    newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.01))
    points(x$domlevels,jitter(x$SR.legumes),col=c, pch=p)
    lines(t(newdata),predict(fit,newdata, se=T)$fit) 
    title(main="Legume richness")
    mtext(3,line=-1, adj=0.1,cex=0.6,
    text=paste('S:', round(cor.test(x$domlevels,x$SR.legumes, method='spearman')$estimate,2), 
     p2star(cor.test(x$domlevels,x$SR.legumes, method='spearman')$p.val)))


i='MELRAM'
x=tmp[as.character(tmp$SpeciesCode)==i,]
    t=rep(1, length(x$SR))
    t[x$vegtype=='W'& !is.na(x$vegtype)]=2
    t[x$vegtype=='G'& !is.na(x$vegtype)]=3
    p=21
    c=c('grey','forestgreen','goldenrod')[t]

       plot(SR.wood~as.numeric(domlevels), x, type='n', xlim=c(1,6))
    fit=lm(SR.wood~domlevels, x)
    newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.01))
    points(x$domlevels,jitter(x$SR.wood),col=c, pch=p)
    lines(t(newdata),predict(fit,newdata, se=T)$fit)
    title(main="Woody richness")
    mtext(3,line=-1, adj=0.1,cex=0.6,
    text=paste('S:', round(cor.test(x$domlevels,x$SR.wood, method='spearman')$estimate,2),
    p2star(cor.test(x$domlevels,x$SR.wood, method='spearman')$p.val)))

    plot(SR.legumes~as.numeric(domlevels), x, type='n', xlim=c(1,6))
    fit=lm(SR.legumes~domlevels, x)
    newdata=data.frame(domlevels=seq(min(x$domlevels, na.rm=T),max(x$domlevels, na.rm=T),0.01))
    points(x$domlevels,jitter(x$SR.legumes),col=c, pch=p)
    lines(t(newdata),predict(fit,newdata, se=T)$fit) 
    title(main="Legume richness")
    mtext(3,line=-1, adj=0.1,cex=0.6,
    text=paste('S:', round(cor.test(x$domlevels,x$SR.legumes, method='spearman')$estimate,2), 
     p2star(cor.test(x$domlevels,x$SR.legumes, method='spearman')$p.val)))

          
mtext(1,text="Abundance class", outer=T, cex=0.8)      
mtext(2,text="Lolium perenne (Alien)", outer=T, at=0.75, line=1.5, cex=0.8, las=0)
mtext(2,text="Melicytus ramiflorus (Native)", outer=T, at=0.25, line=1.5, cex=0.8, las=0)
          
 
############ Negative IMPACTS ############
# Summary table of effects to be exported 
#write.csv(effect.summary[effect.summary$total.effect>0,], file='alleffects.csv')


# species with significant negative impact
negatives=na.omit(as.character(effect.summary$SpeciesCode[effect.summary$effect==-1]))

# Graph
x11()
tmp=databp
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
    text=paste('S:', round(effect.summary$rho[effect.summary$SpeciesCode==i],2),
    p2star(effect.summary$S.pval[effect.summary$SpeciesCode==i])))
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
    text=paste('K:',p2star(effect.summary$k.pval[effect.summary$SpeciesCode==i])))

  }
  }
  
# Graph
  x11()
tmp=databp
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
    text=paste('S:', round(effect.summary$rho[effect.summary$SpeciesCode==i],2),
    p2star(effect.summary$S.pval[effect.summary$SpeciesCode==i])))
    
    mtext(3,line=-1.7, adj=0.1,cex=0.6,
    text=paste('K:',p2star(effect.summary$k.pval[effect.summary$SpeciesCode==i])))
    
    }

## Impact on natives
  x11()
tmp=databp
tmp=tmp[tmp$SpeciesCode %in% negatives,]
#tmp=tmp[tmp$SpeciesCode %in% as.character(effect.summary$SpeciesCode[effect.summary$effect!=-1 & effect.summary$nat.effect==-1]),]

tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(2,5), mar=c(3,2,2,1))
# par(mfrow=c(1,3), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRnat~domlevels, x, type='n', xlim=c(1,7))

  if (effect.summary$nat.effect[effect.summary$SpeciesCode==i]==-1) {
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
        text=paste('S:', round(effect.summary$nat.rho[effect.summary$SpeciesCode==i],2),
                   p2star(effect.summary$nat.S.pval[effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(effect.summary$nat.k.pval[effect.summary$SpeciesCode==i])))
                 
  }
 
## Impact on aliens
  x11()
tmp=databp
tmp=tmp[tmp$SpeciesCode %in% negatives,]
#tmp=tmp[tmp$SpeciesCode %in% as.character(effect.summary$SpeciesCode[effect.summary$ali.effect==-1 & effect.summary$effect>-1]),]

tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(2,5), mar=c(3,2,2,1))
# par(mfrow=c(1,3), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRali~domlevels, x, type='n', xlim=c(1,7))

  if (effect.summary$ali.effect[effect.summary$SpeciesCode==i]==-1) {
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
        text=paste('S:', round(effect.summary$ali.rho[effect.summary$SpeciesCode==i],2),
                   p2star(effect.summary$ali.S.pval[effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(effect.summary$ali.k.pval[effect.summary$SpeciesCode==i])))
                 
  }
 
############ POSITIVE IMPACTS ############

# species with significant positive impact
positives=as.character(effect.summary$SpeciesCode[effect.summary$effect==1])

# Graph
tmp=databp
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
tmp=databp
tmp=tmp[tmp$SpeciesCode %in%positives,]
tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(1,3), mar=c(3,2,2,1))
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
    text=paste('S:', round(effect.summary$rho[effect.summary$SpeciesCode==i],2),
    p2star(effect.summary$S.pval[effect.summary$SpeciesCode==i])))
    
    mtext(3,line=-1.7, adj=0.1,cex=0.6,
    text=paste('K:',p2star(effect.summary$k.pval[effect.summary$SpeciesCode==i])))
    
    }

  ## Impact on natives
  x11()
tmp=databp
tmp=tmp[tmp$SpeciesCode %in% positives,]

tmp=tmp[order(-tmp$ALIEN),]
par(mfrow=c(1,3), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRnat~domlevels, x, type='n', xlim=c(1,7))

  if (effect.summary$nat.effect[effect.summary$SpeciesCode==i]==1) {
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
        text=paste('S:', round(effect.summary$nat.rho[effect.summary$SpeciesCode==i],2),
                   p2star(effect.summary$nat.S.pval[effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(effect.summary$nat.k.pval[effect.summary$SpeciesCode==i])))
                 
  }
 
## Impact on aliens
  x11()
tmp=databp
tmp=tmp[na.omit(tmp$SpeciesCode %in% positives),]
tmp=na.omit(tmp[tmp$SpeciesCode %in% as.character(effect.summary$SpeciesCode[effect.summary$ali.effect==1 & effect.summary$effect<1]),])

tmp=tmp[order(-tmp$ALIEN),]
 par(mfrow=c(2,4), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  plot(SRali~domlevels, x, type='n', xlim=c(1,7))

  if (effect.summary$ali.effect[effect.summary$SpeciesCode==i]==1) {
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
        text=paste('S:', round(effect.summary$ali.rho[effect.summary$SpeciesCode==i],2),
                   p2star(effect.summary$ali.S.pval[effect.summary$SpeciesCode==i])))
  
   mtext(3,line=-1.7, adj=0.1,cex=0.6,
        text=paste('K:',p2star(effect.summary$ali.k.pval[effect.summary$SpeciesCode==i])))
                 
  }