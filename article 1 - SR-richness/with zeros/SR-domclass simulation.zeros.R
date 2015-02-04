# Testing Species richness - species dominance class relationships aginast a simulation
setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R")
######## SIMULATION 1 :  shuffling abundances within communities

# sim1=simul1.class(db=databp,var="SR",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1, file='sim1.Rdata')
# sim1.nat=simul1.class(db=databp,var="SRnat",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.nat, file='sim1.nat.Rdata')
# sim1.ali=simul1.class(db=databp,var="SRali",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.ali, file='sim1.ali.Rdata')
# sim1.grass=simul1.class(db=databp[databp$PlotName%in% grasslands,],var="SR",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.grass, file='sim1.grass.Rdata')
# sim1.wood=simul1.class(db=databp[databp$PlotName%in% woodlands,],var="SR",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.wood, file='sim1.wood.Rdata')
# 
# # sim1.onwood=simul1.class(db=databp,var="SR.wood",min=10, alpha= 0.01, nreps=99)
# # sim1.onlegs=simul1.class(db=databp,var="SR.legumes",min=10, alpha= 0.01, nreps=99)
# # sim1.onlegs.grass=simul1.class(db=databp[databp$PlotName%in% grasslands,],envplot=envplot,zeros=T,var="SR.legumes",min=10, alpha= 0.01, nreps=99)
# # sim1.onlegs.wood=simul1.class(db=databp[databp$PlotName%in% woodlands,],envplot=envplot,zeros=T,var="SR.legumes",min=10, alpha= 0.01, nreps=99)
# 
# sim1.ali.grass=simul1.class(db=databp[databp$PlotName%in% grasslands,],envplot=envplot,zeros=T,var="SRali",min=10, alpha= 0.01, nreps=999); save(sim1.ali.grass, file='sim1.ali.grass.Rdata')
# sim1.ali.wood=simul1.class(db=databp[databp$PlotName%in% woodlands,],envplot=envplot,zeros=T,var="SRali",min=10, alpha= 0.01, nreps=999); save(sim1.ali.wood, file='sim1.ali.wood.Rdata')
# sim1.nat.grass=simul1.class(db=databp[databp$PlotName%in% grasslands,],envplot=envplot,zeros=T,var="SRnati",min=10, alpha= 0.01, nreps=999); save(sim1.nat.grass, file='sim1.nat.grass.Rdata')
# sim1.nat.wood=simul1.class(db=databp[databp$PlotName%in% woodlands,],envplot=envplot,zeros=T,var="SRnat",min=10, alpha= 0.01, nreps=999); save(sim1.nat.wood, file='sim1.nat.wood.Rdata')
# 
### Plotting for results of simulations 1
x11()
# par(mfrow=c(1,3), oma=c(1,1,1,1), mar=c(3,4,1,1), cex=0.9)
names=c("significant","negative.Native","negative.Alien","positive.Native","positive.Alien")
y=data.frame(all=sim1[1:100,names], grass= sim1.grass[1:1000,names],wood=sim1.wood[1:1000,names])

    # significant differences ?
   sig.inf=sapply(1:15,FUN=function(x) {
      p=(sum((y[,x]<y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
      return(p)   
     } )      
                                         
  sig.sup=sapply(1:15,FUN=function(x) {
       p=(sum((y[,x]>y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
      return( p )   
     } )


boxplot(y[,],
        pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5), xaxt='n',yaxt='n',xlim=c(0.5,15.5),
        outline=F, ylim=c(-1.5,round(max(y))+2),col=c('lightgrey',rep(c('white', 'grey30'),2)), border='grey30',whisklty="solid", cex.axis=0.9)
points(1:15, y[1,],  pch=20, cex=1.5)
axis(1, at=sort(c(seq(1,11,5),seq(2.5,12.5,5),seq(4.5,14.5,5))), line=-0.8, labels = rep(c('all','-','+'),3), tcl=-0.3, lty=0, cex.axis=0.8)
axis(1, at=c(3,8,13), line=0.5,lty=0,hadj=0.5,labels = c('All sites','Grasslands', 'Woodlands'),cex.axis=0.8)
axis(1, at=c(0.5,5.5, 10.5,15.5), lwd=0,lwd.ticks=1, labels = NA, tcl=-1)
axis(1, at=sort(c(seq(1.5,11.5,5),seq(3.5,13.5,5))), lwd=0,lwd.ticks=1, labels = NA, tcl=-0.2)
axis(2, las=1, cex.axis=0.8, mgp=c(3,0.7,0), tcl=-0.4 )
title(ylab='Frequency', line=2.5, cex.lab=1)
abline(v=c(5.5, 10.5), lty='dotted')                                                    
box(bty='o')  
sig=apply(cbind(sig.inf, sig.sup),1, min, na.rm=T)
sig=sapply(sig, p2star)
sig[sig=='ns']=''
text(y=y[1,]+2, x=1:15, label=sig, cex.lab=0.7)
mtext(side=1,line=-1, at=1:15, text=y[1,], cex=0.7, col='grey40', font=3)
#        
# mtext(side=3,line=0,adj=0.5,  text="Total species richness", cex=0.9)

legend(x=12, y=14, legend=c('All species', 'Natives','Aliens'), fill=c('lightgrey', 'white','grey30'), cex=0.8, bty='n')

## NAtive richness
names=c("significant","negative.Native","negative.Alien","positive.Native","positive.Alien")
y=data.frame(all=sim1.nat[1:1000,names], grass=sim1.nat.grass[1:1000,names],wood=sim1.nat.wood[1:1000,names])
   # significant differences ?
   sig.inf=sapply(1:15,FUN=function(x) {
      p=(sum((y[,x]<y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
      return(p)   
     } )      
                                         
  sig.sup=sapply(1:15,FUN=function(x) {
       p=(sum((y[,x]>y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
      return( p )   
     } )


boxplot(y[,],
        pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5), xaxt='n',yaxt='n',xlim=c(0.5,15.5),
        outline=F, ylim=c(-1.5,round(max(y, na.rm=T))+2),col=c('lightgrey',rep(c('white', 'grey30'),2)), border='grey30',whisklty="solid", cex.axis=0.9)
points(1:15, y[1,],  pch=20, cex=1.5)
axis(1, at=sort(c(seq(1,11,5),seq(2.5,12.5,5),seq(4.5,14.5,5))), line=-0.8, labels = rep(c('all','-','+'),3), tcl=-0.3, lty=0, cex.axis=0.8)
axis(1, at=c(3,8,13), line=0.5,lty=0,hadj=0.5,labels = c('All sites','Grasslands', 'Woodlands'),cex.axis=0.8)
axis(1, at=c(0.5,5.5, 10.5,15.5), lwd=0,lwd.ticks=1, labels = NA, tcl=-1)
axis(1, at=sort(c(seq(1.5,11.5,5),seq(3.5,13.5,5))), lwd=0,lwd.ticks=1, labels = NA, tcl=-0.2)
axis(2, las=1, cex.axis=0.8, mgp=c(3,0.7,0), tcl=-0.4 )
title(ylab='Frequency', line=2.5, cex.lab=1)
abline(v=c(5.5, 10.5), lty='dotted')                                                    
box(bty='o')  
sig=apply(cbind(sig.inf, sig.sup),1, min, na.rm=T)
sig=sapply(sig, p2star)
sig[sig=='ns']=''
text(y=y[1,]+2, x=1:15, label=sig, cex.lab=0.7)
mtext(side=1,line=-1, at=1:15, text=y[1,], cex=0.7, col='grey40', font=3)

mtext(side=3,line=0,adj=0.5,  text="Native richness", cex=0.9)


# alien richness
names=c("significant","negative.Native","negative.Alien","positive.Native","positive.Alien")
y=data.frame(all=sim1.ali[1:1000,names], grass= sim1.ali.grass[1:1000,names],wood=complete.sim(sim1.ali.wood)[1:1000,names])
  # significant differences ?
   sig.inf=sapply(1:15,FUN=function(x) {
      p=(sum((y[,x]<y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
      return(p)   
     } )      
                                         
  sig.sup=sapply(1:15,FUN=function(x) {
       p=(sum((y[,x]>y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
      return( p )   
     } )


boxplot(y[,],
        pars = list(boxwex = 0.3, staplewex = 0.5, outwex = 0.5), xaxt='n',yaxt='n',xlim=c(0.5,15.5),
        outline=F, ylim=c(-1.5,round(max(y))+2),col=c('lightgrey',rep(c('white', 'grey30'),2)), border='grey30',whisklty="solid", cex.axis=0.9)
points(1:15, y[1,],  pch=20, cex=1.5)
axis(1, at=sort(c(seq(1,11,5),seq(2.5,12.5,5),seq(4.5,14.5,5))), line=-0.8, labels = rep(c('all','-','+'),3), tcl=-0.3, lty=0, cex.axis=0.8)
axis(1, at=c(3,8,13), line=0.5,lty=0,hadj=0.5,labels = c('All sites','Grasslands', 'Woodlands'),cex.axis=0.8)
axis(1, at=c(0.5,5.5, 10.5,15.5), lwd=0,lwd.ticks=1, labels = NA, tcl=-1)
axis(1, at=sort(c(seq(1.5,11.5,5),seq(3.5,13.5,5))), lwd=0,lwd.ticks=1, labels = NA, tcl=-0.2)
axis(2, las=1, cex.axis=0.8, mgp=c(3,0.7,0), tcl=-0.4 )
title(ylab='Frequency', line=2.5, cex.lab=1)
abline(v=c(5.5, 10.5), lty='dotted')                                                    
box(bty='o')  
sig=apply(cbind(sig.inf, sig.sup),1, min, na.rm=T)
sig=sapply(sig, p2star)
sig[sig=='ns']=''
text(y=y[1,]+2, x=1:15, label=sig, cex.lab=0.7)
mtext(side=1,line=-1, at=1:15, text=y[1,], cex=0.7, col='grey40', font=3)
mtext(side=3,line=0,adj=0.5,  text="Alien richness", cex=0.9)
        
        
              
################################################################################################################  

######## SIMULATION 1 : rho values only !!
       
sim1.rhos.zeros=simul1.class.rhos(db=databp,var="SR",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.rhos, file='sim1.rhos.zeros.Rdata')
sim1.rhos.zeros.grass=simul1.class.rhos(db=databp[databp$PlotName%in% grasslands,],var="SR",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.rhos.zeros.grass, file='sim1.rhos.zeros.grass.Rdata')
sim1.rhos.zeros.wood=simul1.class.rhos(db=databp[databp$PlotName%in% woodlands,],var="SR",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.rhos.zeros.wood, file='sim1.rhos.zeros.wood.Rdata')

sim1.rhos.zeros.nat=simul1.class.rhos(db=databp,var="SRnat",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.rhos.zeros.nat, file='sim1.rhos.zeros.nat.Rdata')
sim1.rhos.zeros.ali=simul1.class.rhos(db=databp,var="SRali",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.rhos.zeros.ali, file='sim1.rhos.zeros.ali.Rdata')

sim1.rhos.zeros.ali.grass=simul1.class.rhos(db=databp[databp$PlotName%in% grasslands,],var="SRali",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.rhos.zeros.ali.grass, file='sim1.rhos.zeros.ali.grass.Rdata')
sim1.rhos.zeros.ali.wood=simul1.class.rhos(db=databp[databp$PlotName%in% woodlands,],var="SRali",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.rhos.zeros.ali.wood, file='sim1.rhos.zeros.ali.wood.Rdata')
sim1.rhos.zeros.nat.grass=simul1.class.rhos(db=databp[databp$PlotName%in% grasslands,],var="SRnat",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.rhos.zeros.nat.grass, file='sim1.rhos.zeros.nat.grass.Rdata')
sim1.rhos.zeros.nat.wood=simul1.class.rhos(db=databp[databp$PlotName%in% woodlands,],var="SRnat",envplot=envplot,zeros=T,min=10, alpha= 0.01, nreps=999); save(sim1.rhos.zeros.nat.wood, file='sim1.rhos.zeros.nat.wood.Rdata')

shuff.all=shuff.classSR(db = databp,species = species, envplot=envplot,zeros=T,min.occur = 10, min.class = 3, alpha = 0.01, nreps = 999) ; save(shuff.all, file='shuff.all.Rdata')
shuff.grass=shuff.classSR(db = databp[databp$PlotName%in% grasslands,],species = species, min.occur = 10, envplot=envplot,zeros=T,min.class = 3, alpha = 0.01, nreps = 999) ; save(shuff.grass, file='shuff.grass.Rdata')
shuff.wood=shuff.classSR(db = databp[databp$PlotName%in% woodlands,],species = species, min.occur = 10, envplot=envplot,zeros=T,min.class = 3, alpha = 0.01, nreps = 999) ; save(shuff.wood, file='shuff.wood.Rdata')



### Boxplots and ttests with null values
x11()
par(mfrow=c(3,3), mar=c(2,1,1,0),oma=c(1,8,1,1), las=1, cex=0.9)

results.simul1.rhos(x=sim1.rhos)
axis(2, cex.axis=0.7, las=1,mgp=c(3, 0.7, 0))
mtext(side=2, text=expression("Spearman's "*rho), outer=F, las=0, line=1.7, cex=0.7)

results.simul1.rhos(x=sim1.rhos.zeros.nat)
results.simul1.rhos(x=sim1.rhos.zeros.ali)

results.simul1.rhos(x=sim1.rhos.zeros.grass)
axis(2, cex.axis=0.7, las=1,mgp=c(3, 0.7, 0))
mtext(side=2, text=expression("Spearman's "*rho), outer=F, las=0, line=1.7, cex=0.7)

results.simul1.rhos(x=sim1.rhos.zeros.nat.grass)
results.simul1.rhos(x=sim1.rhos.zeros.ali.grass)

results.simul1.rhos(x=sim1.rhos.zeros.wood)
axis(2, cex.axis=0.7, las=1,mgp=c(3, 0.7, 0))
mtext(side=2, text=expression("Spearman's "*rho), outer=F, las=0, line=1.7, cex=0.7)

results.simul1.rhos(x=sim1.rhos.zeros.nat.wood)
results.simul1.rhos(x=sim1.rhos.zeros.ali.wood)

mtext(side=3, at=c(0.18,0.52,0.85), text=c('Total richness', 'Native richness', 'Alien richness'),
      outer=T, las=0, line=-0.5, cex=0.8)

# mtext(side=2, at=c(0.18,0.52,0.85), text=c('Grasslands', 'Woodlands','All plots'),
#       outer=T, las=1, line=2.5, adj=1)
mtext(side=2, at=c(0.28,0.62,0.95), text=c( 'Woodlands','Grasslands','All plots'),
      outer=T, las=1, line=2, adj=1, font=3, cex=0.8)


#### autre orga
x11()
par(mfrow=c(3,3), mar=c(2,1,0,0),oma=c(1,8,3,1), las=1, cex=0.9)

results.simul1.rhos(x=sim1.rhos)
axis(2, cex.axis=0.7, las=1,mgp=c(3, 0.7, 0))
mtext(side=2, text=expression("Spearman's "*rho), outer=F, las=0, line=2, cex=0.7)
results.simul1.rhos(x=sim1.rhos.zeros.grass)
results.simul1.rhos(x=sim1.rhos.zeros.wood)

results.simul1.rhos(x=sim1.rhos.zeros.nat)

axis(2, cex.axis=0.7, las=1,mgp=c(3, 0.7, 0))
mtext(side=2, text=expression("Spearman's "*rho), outer=F, las=0, line=1.7, cex=0.7)
results.simul1.rhos(x=sim1.rhos.zeros.nat.grass)
results.simul1.rhos(x=sim1.rhos.zeros.nat.wood)

results.simul1.rhos(x=sim1.rhos.zeros.ali)
axis(2, cex.axis=0.7, las=1,mgp=c(3, 0.7, 0))
mtext(side=2, text=expression("Spearman's "*rho), outer=F, las=0, line=1.7, cex=0.7)
results.simul1.rhos(x=sim1.rhos.zeros.ali.grass)
results.simul1.rhos(x=sim1.rhos.zeros.ali.wood)

mtext(side=3, at=c(0.18,0.52,0.85), text=c('All plots', 'Grasslands','Woodlands'),
      outer=T, las=0, line=0.5, cex=0.9)

# mtext(side=2, at=c(0.21,0.54,0.88), text=c('c)  Alien\n     richness', 'b)  Native\n      richness', 'a)  Total\n      richness'),
#       outer=T, las=1, line=6.5, adj=0, font=0, cex=0.9)

# mtext(side=2, at=c(0.27,0.60,0.94), text=c('c)', 'b)', 'a)'),
#       outer=T, las=1, line=4, adj=0, font=0, cex=0.9)

mtext(side=2, at=c(0.21,0.54,0.88), text=c('Alien\nrichness', 'Native\nrichness', 'Total\nrichness'),
      outer=T, las=1, line=4, adj=0.5, font=0, cex=0.9)
mtext(side=2, at=c(0.15,0.48,0.82), text=c('(c)', '(b)', '(a)'),
      outer=T, las=1, line=4, adj=0.54, font=0, cex=0.8)


# mtext(side=2, at=c(0.21,0.54,0.88), text=c('Alien\nrichness\n(c)', 'Native\nrichness\n(b)', 'Total\nrichness\n(a)'),
#       outer=T, las=1, line=4, adj=0.5, font=0, cex=0.9)


## Testing the frequency of combinations of - and + effects
 
# sim1.combi=simul1.combinations(db=databp,min=10, alpha= 0.01, nreps=9999); save(sim1.combi, file='sim1.combi.Rdata')
# sim1.combi.grass=simul1.combinations(db=databp[databp$PlotName%in% grasslands,],min=10, alpha= 0.01, nreps=999); save(sim1.combi.grass, file='sim1.combi.grass.Rdata')
# sim1.combi.wood=simul1.combinations(db=databp[databp$PlotName%in% woodlands,],min=10, alpha= 0.01, nreps=999); save(sim1.combi.wood, file='sim1.combi.wood.Rdata')

#filtered
filtered.freq.combi=results.sim1.combination(sim=sim1.combi, filter=T, by.status=T)

x11()
par(mfrow=c(1,3), mar=c(4, 1,1,1), oma=c(1,4,1,1))
freq.combi=results.sim1.combination(sim=sim1.combi, filter=F, by.status=T)
freq.combi.grass=results.sim1.combination(sim=sim1.combi.grass, filter=F,by.status=T)
freq.combi.wood=results.sim1.combination(sim=sim1.combi.wood, filter=F,by.status=T)
mtext(side=2,text="Number of species",outer=T, las=0, cex=0.8, line=2)
