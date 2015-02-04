### Figures for article abundance richness relationships
load(file="saved Rdata/article 1/JEcol resubmission.Rdata")

############### FUNCTION: illustration for a focal species species #######################################################
plot.illustration=function(i='LOLPER',var="SR", db=envplot, X=comm, zeros=F, 
                           c="grey40",xaxis=T, yaxis=T, vegtype=F) {
  x=data.frame(A=X[match(as.character(db$PLOTID),rownames(X)),i],var=na.omit(as.vector(db[,var])))
  if (zeros)  x$A[x$A==0]=0
  if (!zeros)  x$A[x$A==0]=NA
  
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
  plot(var ~ A, data=x, type='n', xlim=c(0,6), ylim=c(0,70), xaxt='n', yaxt='n', bty='l')
  #     if (xaxis) axis(1, at=1:6 , label=1:6, tcl=0.3, mgp=c(3,0.1,0),cex.axis=0.9)
  #     if (yaxis) axis(2, las=1, cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
  if (xaxis==F) axis(1, at=1:6 , label=NA, tcl=0.3, mgp=c(3,0.1,0),cex.axis=1)
  if (yaxis==F) axis(2, las=1, label=NA,cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
  
  if (xaxis==T) axis(1, at=1:6 , label=1:6, tcl=0.3, mgp=c(3,0.1,0),cex.axis=1)
  if (yaxis==T) axis(2, las=1, cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
  
  fit=lm(var ~ A, data=x)
  newdata=data.frame(domlevels=seq(min(x$A, na.rm=T),max(x$A, na.rm=T),0.01))
  points(x$A,jitter(x$var),col=c, pch=p,bg=bg)
  
  ## spearman rank test correlation 
  co=cor.test(x$A, x$var, method="spearman")
  mtext(3,line=-1.1, adj=1,cex=0.7,
        text=substitute(italic(rho)*' = '*x, list(x=paste(round(co$estim,2),p2star(co$p.val), sep=''))))
  if(co$p.val<0.01) lines(na.omit(x$A),predict(fit))
}

############### Figure 1 ###################################
par(mfcol=c(3,4), oma=c(3,10,4,1), mar=c(1,1,0,0), las=1)
plot.illustration(i='CRIMUR',var="SR", xaxis=F, yaxis=T)
plot.illustration(i='CRIMUR',var="SRnat", xaxis=F, yaxis=T)
plot.illustration(i='CRIMUR',var="SRali", xaxis=T, yaxis=T)

plot.illustration(i='RUMACE',var="SR", xaxis=F, yaxis=F)
plot.illustration(i='RUMACE',var="SRnat", xaxis=F, yaxis=F)
plot.illustration(i='RUMACE',var="SRali", xaxis=T, yaxis=F)

plot.illustration(i='GRILIT',var="SR", xaxis=F, yaxis=F)
plot.illustration(i='GRILIT',var="SRnat", xaxis=F, yaxis=F)
plot.illustration(i='GRILIT',var="SRali", xaxis=T, yaxis=F)

plot.illustration(i='COPROT',var="SR", xaxis=F, yaxis=F)
plot.illustration(i='COPROT',var="SRnat", xaxis=F, yaxis=F)
plot.illustration(i='COPROT',var="SRali", xaxis=T, yaxis=F)

mtext(1,text="Abundance class", outer=T, line=0.5,cex=0.9)      
mtext(3,text=expression(italic("a) C. murinum\n   (Alien grass)")), outer=T, at=0.12, line=1, cex=0.8, las=1, adj=0.5)
mtext(3,text=expression(italic("b) R. acetosella\n   (Alien herb)")), outer=T, at=0.37, line=1, cex=0.8, las=1, adj=0.5)
mtext(3,text=expression(italic("c) G. littoralis\n   (Native tree)")), outer=T, at=0.62, line=1, cex=0.8, las=1, adj=0.5)
mtext(3,text=expression(italic("d) C. rotundifolia\n   (Native shrub)")), outer=T, at=0.87, line=1, cex=0.8, las=1, adj=0.5)

mtext(2,text="Species richness", outer=T, at=0.52, line=8, cex=0.9, las=0)

mtext(2,text="Alien\nrichness", outer=T, at=0.18, line=2, cex=0.9, las=1)
mtext(2,text="Native\nrichness", outer=T, at=0.51, line=2, cex=0.9, las=1)
mtext(2,text="Total\nrichness", outer=T, at=0.85, line=2, cex=0.9, las=1)


# FUNCTION: PLOT RHO DISTRIBUTIONS (null model 1) #######################################################
 plot.simul.rhos=function(x, ylab=T){
  plot(c(0.5,3.5),c(min(x[,1], na.rm=T),max(x[,1], na.rm=T)), type='n', ann=F, axes=F,ylim=c(-0.8,0.8))
  
  abline(h=0)
  quant=cbind(quantile(x[,1], probs=c(0.0275,0.9725), na.rm=T),quantile(x[,3], probs=c(0.0275,0.9725), na.rm=T),quantile(x[,5], probs=c(0.0275,0.9725), na.rm=T))
  arrows(1:3, quant[1,],1:3, quant[2,], col='black',angle = 90,code =3, length=0.05)  

  tmp=boxplot(x[,c(1,3,5)], at=1:3, ann=F,axes=F,col=c('lightgrey', 'white', 'grey40'),
           pars = list(boxwex = 0.3, staplewex = 0, outwex = 0),
         whisklty='blank',  staplelty='blank', boxlty='solid',outline=F,add=T,
          names=NULL,xaxt='n', cex.axis=0.8, las=1)
#  points(1:3,apply(x[,c(1,3,5)],2,FUN=median),  pch=21, cex=1.5, bg=c('black','grey80','grey40') )

  out=rbind(cbind(x$rho[which(x$rho<quant[1,1] | x$rho>quant[2,1])],1),
            cbind(x$NR.rho[which(x$NR.rho<quant[1,2] | x$NR.rho>quant[2,2])],2) ,
            cbind(x$AR.rho[which(x$AR.rho<quant[1,3] | x$AR.rho>quant[2,3])],3))
  points(out[,2],out[,1],pch=".", cex=1)

#   axis(1,at=1:3, label=c('Total\nrichness','Native\nrichness', 'Alien\nrichness'), cex.axis=0.8, tcl=-0.3, mgp=c(3,0.8,0))

   axis(2, at=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75),label= c(NA,-0.5,NA,0,NA,0.5,NA),
         cex.axis=0.7, tcl=0.18,las=1,mgp=c(3, 0.3, 0))

#   box(bty='o')

# TEST difference from null model mean
y= quantile(x$rho, 0.75, na.rm=T) +0.1
t=wilcox.test(x$rho,x[,2], paired=T)
text(x=1.2, y=y, label=(p2star(t$p.val)), cex=0.9)
 
y= quantile(x$NR.rho, 0.75, na.rm=T) +0.1
t=wilcox.test(x$NR.rho,x[,4], paired=T)
text(x=2.2, y=y, label=(p2star(t$p.val)), cex=0.9)
 
y= quantile(x$AR.rho, 0.75, na.rm=T) +0.1
t=wilcox.test(x$AR.rho,x[,6], paired=T)
text(x=3.2, y=y, label=(p2star(t$p.val)), cex=0.9) 
 
if (ylab) mtext(side=2, text=expression("Spearman's "*italic(rho)), outer=F, las=0, line=1.2, cex=0.75)

}

####### FIG 2 : rhos in ALL PLOTS ##########
x=data.frame( rho= sim1.rhos[,1], mean.null= rowMeans(sim1.rhos, na.rm=T),
              NR.rho= sim1.rhos.nat[,1], NR.mean.null= rowMeans(sim1.rhos.nat, na.rm=T),
              AR.rho= sim1.rhos.ali[,1], AR.mean.null= rowMeans(sim1.rhos.ali, na.rm=T))
x$ALIEN= species[rownames(sim1.rhos), "ALIEN"]
x$WEED= species[rownames(sim1.rhos), "WEEDOC"]

quartz()

par(mfcol=c(3,1), mar=c(0,0,0,0),oma=c(2,7,2,0), las=1, cex=0.9)

plot.simul.rhos(x)
# axis(3,at=1:3, label=c('Total\nrichness','Native\nrichness', 'Alien\nrichness'), 
#      cex.axis=0.75,lwd=0, tcl=-0.3, mgp=c(3,-0.5,0))
axis(3,at=1:3, label=c('Total\nrichness','Native\nrichness', 'Alien\nrichness'), 
     cex.axis=0.8,lwd=0, tcl=-0.3, mgp=c(3,-0.2,0), font=3)

mtext(text="a)", side=3, line=-1.4, adj=0.05, cex=0.8)

plot.simul.rhos(x[x$ALIEN==0,])
mtext(text="b)", side=3, line=-1.4, adj=0.05, cex=0.8)
plot.simul.rhos(x[x$ALIEN==1,]) 
# plot.simul.rhos(x[x$ALIEN==1 & x$WEED==0 ,])
# plot.simul.rhos(x[x$ALIEN==1 & x$WEED==1 ,])
mtext(text="c)", side=3, line=-1.4, adj=0.05, cex=0.8)    
#      mtext(side=2, at=c(0.15,0.48,0.82), text=c('c)', 'b)', 'a)'),
#       outer=T, las=1, line=4, adj=0.54, font=0, cex=0.8)
#       mtext(side=2, at=c(0.21,0.54,0.88), text=c('Alien focal species', 'Native focal species','All species'),
#       outer=T, las=0, line=3, cex=0.8, las=1)
mtext(side=2, at=c(0.18,0.52,0.84), text=c('Alien\nfocal species', 'Native\nfocal species','All\nfocal species'),
      outer=T, las=0, line=4.6, cex=0.8, adj=0.5,las=1)

####### FIG 3 : Rhos in woodland and grasslands  ##########

quartz()
par(mfcol=c(3,2), mar=c(0,3,0,0),oma=c(0,5,4,0), las=1, cex=0.9)
     
# Grasslands
x=data.frame( rho= sim1.rhos.grass[,1], mean.null= rowMeans(sim1.rhos.grass, na.rm=T),
              NR.rho= sim1.rhos.nat.grass[,1], NR.mean.null= rowMeans(sim1.rhos.nat.grass, na.rm=T),
              AR.rho= sim1.rhos.ali.grass[,1], AR.mean.null= rowMeans(sim1.rhos.ali.grass, na.rm=T))
x$ALIEN= species[rownames(sim1.rhos.grass), "ALIEN"]
x$WEED= species[rownames(sim1.rhos.grass), "WEED"]

plot.simul.rhos(x)
axis(3,at=1:3, label=c('Total\nrichness','Native\nrichness', 'Alien\nrichness'), 
     cex.axis=0.8,lwd=0, tcl=-0.3, mgp=c(3,-0.2,0), font=3)
# mtext(text="a)", side=3, line=-1, adj=-0.1, cex=0.8)
mtext(text="a)", side=3, line=-1.4, adj=0.05, cex=0.8)
 plot.simul.rhos(x[x$ALIEN==0,])
 mtext(text="b)", side=3, line=-1.4, adj=0.05, cex=0.8)
plot.simul.rhos(x[x$ALIEN==1,])
 mtext(text="c)", side=3, line=-1.4, adj=0.05, cex=0.8)
# plot.simul.rhos(x[x$ALIEN==1 & x$WEED==0 ,])
# plot.simul.rhos(x[x$ALIEN==1 & x$WEED==1 ,])
 
# Woodlands
x=data.frame( rho= sim1.rhos.wood[,1], mean.null= rowMeans(sim1.rhos.wood, na.rm=T),
              NR.rho= sim1.rhos.nat.wood[,1], NR.mean.null= rowMeans(sim1.rhos.nat.wood, na.rm=T),
              AR.rho= sim1.rhos.ali.wood[,1], AR.mean.null= rowMeans(sim1.rhos.ali.wood, na.rm=T))
x$ALIEN= species[rownames(sim1.rhos.wood), "ALIEN"]
x$WEED= species[rownames(sim1.rhos.wood), "WEED"]

plot.simul.rhos(x, ylab=F)
axis(3,at=1:3, label=c('Total\nrichness','Native\nrichness', 'Alien\nrichness'), 
     cex.axis=0.8,lwd=0, tcl=-0.3, mgp=c(3,-0.2,0), font=3)
 mtext(text="d)", side=3, line=-1.4, adj=0.05, cex=0.8)

 plot.simul.rhos(x[x$ALIEN==0,], ylab=F)
  mtext(text="e)", side=3, line=-1.4, adj=0.05, cex=0.8)
 plot.simul.rhos(x[x$ALIEN==1,], ylab=F)
  mtext(text="f)", side=3, line=-1.4, adj=0.05, cex=0.8)
# plot.simul.rhos(x[x$ALIEN==1 & x$WEED==0 ,])
# plot.simul.rhos(x[x$ALIEN==1 & x$WEED==1 ,])
     
# mtext(side=2, at=c(0.15,0.48,0.82), text=c('(c)', '(b)', '(a)'),
# outer=T, las=1, line=3, adj=0.54, font=0, cex=0.8)
mtext(side=2, at=c(0.18,0.52,0.84), text=c('Alien\nfocal species', 'Native\nfocal species','All\nfocal species'),
outer=T, las=0, line=2, cex=0.8, adj=0.5,las=1)

mtext(side=3, at=c(0.3, 0.80), text=c('Grasslands','Woodlands'),
      outer=T, las=1, line=2.5, adj=0.5, font=0, cex=1)

     


# FUNCTION: PLOT SIGNIFICANT SPECIES (null model 2) #######################################################
plot.sim2=function(y, space=0.9){  
  x=c(0.9,1.1,1.9,2.1, 2.9, 3.1)
  
    # significant differences ?
       sig.inf=sapply(1:6,FUN=function(x) {
          p=(sum((y[,x]<y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
          return(p)   
         } )      
                                             
      sig.sup=sapply(1:6,FUN=function(x) {
           p=(sum((y[,x]>y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
          return( p )   
         } )  
  
    # SES
    ses=sapply(1:6,FUN=function(x) {
          ses= (y[1,x] - mean(y[,x], na.rm=T)) /sd(y[,x], na.rm=T)
          return(ses)   
         } )  
    ses[is.nan(ses)] =NA  
    f=c(y[1,1],-y[1,2],y[1,3],-y[1,4],y[1,5],-y[1,6])
    nf=c(mean(y[,1]),-mean(y[,2]),mean(y[,3]),-mean(y[,4]),mean(y[,5]),-mean(y[,6]))
   nf.sd =c(sd(y[,1]),-sd(y[,2]),sd(y[,3]),-sd(y[,4]),sd(y[,5]),-sd(y[,6]))
  # initiate plot
    plot(x, f, type='n',  col=c('grey','black'), ann=F, axes=F, xlim=c(0.5,3.5),
         ylim=ylim)
    
  # observed values  
  par(new=T)    
  polygon(x=c(0.80,0.80,1,1), c(f[1],0,0,f[1]), col='lightgrey') 
  polygon(x=c(1,1,1.20,1.20),c(0, f[2],f[2],0),  col='lightgrey')

  polygon(x=c(1.80,1.80,2,2), c(f[3],0,0,f[3]), col='white')
  polygon(x=c(2,2,2.20,2.20), c(0, f[4],f[4],0), col='white')

  polygon(x=c(2.80,2.80,3,3), c(f[5],0,0,f[5]), col='grey40')
  polygon(x=c(3,3,3.20,3.20), c(0, f[6],f[6],0), col='grey40')

#     sig=apply(cbind(sig.inf, sig.sup),1, min, na.rm=T)
    sig=sapply(sig.sup, p2star)
    sig[sig=='ns']=''
 
  abline(h=0)
  ##STARS :
       stars=sapply(1:6, FUN=function(i)  max(abs(f[i]),abs(nf[i]+nf.sd[i]), na.rm=T))
       stars=stars*c(1,-1)
      if(sum(stars>0, na.rm=T)>0) text(y=stars[stars>0]+2, x=x[stars>0], label=sig[stars>0], cex.lab=1)
      if(sum(stars<0, na.rm=T)>0) text(y=stars[stars<0]-2, x=x[stars<0], label=sig[stars<0], cex.lab=1)
      
       axis(2, at=c(-10,-5,0,5,10), label=c(10,5,0,5,10), las=1, cex.axis=1, mgp=c(3,0.2,0), tcl=0.2 )
      axis(2, at=seq(-10,10,1), label=rep("",21), las=1, cex.axis=1, mgp=c(3,0.2,0),lwd=0, lwd.ticks=1, tcl=0.1 )

}

####### FIG 4 : Significant correlations in native and alien focal species  ###############################
 # graph by alien vs natives
 par(mfrow=c(2,3), oma=c(1,10,5,1), mar=c(0,1,0,0),  las=1)
 ylim=c(-13,12)
 
 # all plots
 names=c("Natives.posSR","Natives.negSR","Natives.posSRnat" ,"Natives.negSRnat","Natives.posSRali", "Natives.negSRali")
 plot.sim2(y=shuff.all[1:1000,names])
 axis(3, at=c(1,2,3), line=0,adj=0.5,padj=1, font=3,
      labels = c('Total\nrichness','Native\nrichness','Alien\nrichness'), tcl=-0.3, lty=0, cex.axis=1)
 mtext(2,at=c(5,-5), text=c('+', '-'), outer=F, adj=0.5,line=1.7, cex=1.5, las=1)
 plot.sim2(y=shuff.grass[1:1000,names])
 axis(3, at=c(1,2,3), line=0,adj=0.5,padj=1, font=3,
      labels = c('Total\nrichness','Native\nrichness','Alien\nrichness'), tcl=-0.3, lty=0, cex.axis=1)
 
 plot.sim2(y=shuff.wood[1:1000,names])
 axis(3, at=c(1,2,3), line=0,adj=0.5,padj=1, font=3,
      labels = c('Total\nrichness','Native\nrichness','Alien\nrichness'), tcl=-0.3, lty=0, cex.axis=1)
 
 
 names=c("Aliens.posSR","Aliens.negSR","Aliens.posSRnat" ,"Aliens.negSRnat","Aliens.posSRali", "Aliens.negSRali")
 plot.sim2(y=shuff.all[1:1000,names])
 mtext(2,at=c(5,-5), text=c('+', '-'), outer=F, adj=0.5,line=1.7, cex=1.5, las=1)
 plot.sim2(y=shuff.grass[1:1000,names])
 plot.sim2(y=shuff.wood[1:1000,names])
  
 
 mtext(side=3, at=c(0.18,0.51,0.84),text=c('a) All plots', 'b) Grasslands','c) Woodlands'), 
       outer=T, las=1, line=2, adj=0.5, font=0, cex=0.9)
 mtext(side=2, at=c(0.75,0.25), text=c('Native\nfocal species', 'Alien\nfocal species'),
       outer=T, las=1, line=3.5, cex=0.8, adj=0.5)    
 mtext(2, at=c(0.5),text='Number of significant correlations', outer=T, line=7.5,cex=0.8, las=0)
 

######## graphical abstract #######
plot.illustration=function(i='LOLPER',var="SR", db=envplot, X=comm, zeros=F, 
                           c="grey40",xaxis=T, yaxis=T, vegtype=F) {
  x=data.frame(A=X[match(as.character(db$PLOTID),rownames(X)),i],var=na.omit(as.vector(db[,var])))
  if (zeros)  x$A[x$A==0]=0
  if (!zeros)  x$A[x$A==0]=NA
  
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
  plot(var ~ A, data=x, type='n', xlim=c(0,6), ylim=c(0,50), xaxt='n', yaxt='n', bty='l')
  #     if (xaxis) axis(1, at=1:6 , label=1:6, tcl=0.3, mgp=c(3,0.1,0),cex.axis=0.9)
  #     if (yaxis) axis(2, las=1, cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
  if (xaxis==F) axis(1, at=1:6 , label=NA, tcl=0.3, mgp=c(3,0.1,0),cex.axis=1)
  if (yaxis==F) axis(2, las=1, label=NA,cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
  
  if (xaxis==T) axis(1, at=1:6 , label=1:6, tcl=0.3, mgp=c(3,0.1,0),cex.axis=1)
  if (yaxis==T) axis(2, las=1, cex.axis=1, mgp=c(3,0.3,0), tcl=0.3)
  
  fit=lm(var ~ A, data=x)
  newdata=data.frame(domlevels=seq(min(x$A, na.rm=T),max(x$A, na.rm=T),0.01))
  points(x$A,jitter(x$var),col=c, pch=p,bg=bg)
  
  ## spearman rank test correlation 
  co=cor.test(x$A, x$var, method="spearman")
#   mtext(3,line=-1.1, adj=1,cex=0.7,
#         text=substitute(x, list(x=paste(round(co$estim,2),p2star(co$p.val), sep=''))))
  if(co$p.val<0.01) lines(na.omit(x$A),predict(fit))
}

par(mfcol=c(2,2), oma=c(2,5,2,1), mar=c(1,1,0,0), las=1)
plot.illustration(i='CRIMUR',var="SRnat", xaxis=F, yaxis=T)
plot.illustration(i='CRIMUR',var="SRali", xaxis=T, yaxis=T)

plot.illustration(i='COPROT',var="SRnat", xaxis=F, yaxis=F)
plot.illustration(i='COPROT',var="SRali", xaxis=T, yaxis=F)

mtext(1,text="Species abundance", outer=T, line=0.5,cex=0.9)      
# mtext(3,text=expression("a) Alien grass" *italic("\n (C. murinum)")), outer=T, at=0.25, line=1, cex=1, las=1, adj=0.5)
# mtext(3,text=expression("b) Native shrub "*italic("\n(C. rotundifolia)")), outer=T, at=0.80, line=1, cex=1, las=1, adj=0.5)
mtext(3,text=expression("Alien grass"), outer=T, at=0.25, line=1, cex=1, las=1, padj=0.5, adj=0.5)
mtext(3,text=expression("Native shrub"), outer=T, at=0.80, line=1, cex=1, las=1,padj=0.5, adj=0.5)


mtext(2,text="Alien\nrichness", outer=T, at=0.3, line=1, cex=1, las=1)
mtext(2,text="Native\nrichness", outer=T, at=0.8, line=1, cex=1, las=1)



############## graphical abstract 2 #### 

# FUNCTION: PLOT SIGNIFICANT SPECIES (null model 2) #######################################################
plot.sim2=function(y, space=0.9){  
  x=c(0.9,1.1,1.9,2.1, 2.9, 3.1)
  
  # significant differences ?
  sig.inf=sapply(1:6,FUN=function(x) {
    p=(sum((y[,x]<y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
    return(p)   
  } )      
  
  sig.sup=sapply(1:6,FUN=function(x) {
    p=(sum((y[,x]>y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
    return( p )   
  } )  
  
  # SES
  ses=sapply(1:6,FUN=function(x) {
    ses= (y[1,x] - mean(y[,x], na.rm=T)) /sd(y[,x], na.rm=T)
    return(ses)   
  } )  
  ses[is.nan(ses)] =NA  
  f=c(y[1,1],-y[1,2],y[1,3],-y[1,4],y[1,5],-y[1,6])
  nf=c(mean(y[,1]),-mean(y[,2]),mean(y[,3]),-mean(y[,4]),mean(y[,5]),-mean(y[,6]))
  nf.sd =c(sd(y[,1]),-sd(y[,2]),sd(y[,3]),-sd(y[,4]),sd(y[,5]),-sd(y[,6]))
  # initiate plot
  plot(x, f, type='n',  col=c('grey','black'), ann=F, axes=F, xlim=c(0.5,3.5),
       ylim=ylim)
  
  # observed values  
  par(new=T)    
  polygon(x=c(0.80,0.80,1,1), c(f[1],0,0,f[1]), col='lightgrey') 
  polygon(x=c(1,1,1.20,1.20),c(0, f[2],f[2],0),  col='lightgrey')
  
  polygon(x=c(1.80,1.80,2,2), c(f[3],0,0,f[3]), col='white')
  polygon(x=c(2,2,2.20,2.20), c(0, f[4],f[4],0), col='white')
  
  polygon(x=c(2.80,2.80,3,3), c(f[5],0,0,f[5]), col='grey40')
  polygon(x=c(3,3,3.20,3.20), c(0, f[6],f[6],0), col='grey40')
  
  #     sig=apply(cbind(sig.inf, sig.sup),1, min, na.rm=T)
  sig=sapply(sig.sup, p2star)
  sig[sig=='ns']=''
  
  abline(h=0)
  ##STARS :
  stars=sapply(1:6, FUN=function(i)  max(abs(f[i]),abs(nf[i]+nf.sd[i]), na.rm=T))
  stars=stars*c(1,-1)
  if(sum(stars>0, na.rm=T)>0) text(y=stars[stars>0]+2, x=x[stars>0], label=sig[stars>0], cex.lab=1)
  if(sum(stars<0, na.rm=T)>0) text(y=stars[stars<0]-2, x=x[stars<0], label=sig[stars<0], cex.lab=1)
  
  axis(2, at=c(-10,-5,0,5,10), label=c(10,5,0,5,10), las=1, cex.axis=1, mgp=c(3,0.2,0), tcl=0.2 )
  axis(2, at=seq(-10,10,1), label=rep("",21), las=1, cex.axis=1, mgp=c(3,0.2,0),lwd=0, lwd.ticks=1, tcl=0.1 )
  
}

####### FIG 4 : Significant correlations in native and alien focal species  ###############################
# graph by alien vs natives
par(mfrow=c(1,1), oma=c(0,0,0,0),mar=c(0,8,2,1),  las=1)
ylim=c(-13,12)

# all plots
names=c("Aliens.posSR","Aliens.negSR","Aliens.posSRnat" ,"Aliens.negSRnat","Aliens.posSRali", "Aliens.negSRali")
plot.sim2(y=shuff.all[1:1000,names])
mtext(2,at=c(7,-7), text=c('+', '-'), outer=F, adj=0.5,line=1.7, cex=1.5, las=1)

mtext(2, at=c(0.5),text='Correlations with\nalien species\nabundance', outer=F, line=1.5,cex=0.9, las=1)
axis(3, at=c(1,2,3), line=0,adj=0.5,padj=1, font=3,
     labels = c('Total\nrichness','Native\nrichness','Alien\nrichness'), tcl=-0.3, lty=0, cex.axis=1)

