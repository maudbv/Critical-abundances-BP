# Pairwise abundance classes tests

# source('script/data/import fede data.R')
# source('script/functions/function paired-class SR.R')

# Calculate t-tests between pairs of abundance classes :
pclassSR=pairedclass.test()
pclassSRnat=pairedclass.test(var="SRnat")
pclassSRali=pairedclass.test(var="SRali")

pclassSR.grass=pairedclass.test(db=databp[databp$PlotName %in% grasslands,])
pclassSRnat.grass=pairedclass.test(db=databp[databp$PlotName %in% grasslands,],var="SRnat")
pclassSRali.grass=pairedclass.test(db=databp[databp$PlotName %in% grasslands,],var="SRali")

pclassSR.wood=pairedclass.test(db=databp[databp$PlotName %in% woodlands,])
pclassSRnat.wood=pairedclass.test(db=databp[databp$PlotName %in% woodlands,],var="SRnat")
pclassSRali.wood=pairedclass.test(db=databp[databp$PlotName %in% woodlands,],var="SRali")

## diovided by low and high grasslands

db=databp[databp$PlotName %in% grasslands & databp$DEM_10 < 400,]
length(unique(db$PlotName))

pclassSR.low=pairedclass.test(db=db)
pclassSRnat.low=pairedclass.test(db=db,var="SRnat")
pclassSRali.low=pairedclass.test(db=db,var="SRali")

db=databp[databp$PlotName %in% grasslands & databp$DEM_10 >= 400,]
length(unique(db$PlotName))
pclassSR.high=pairedclass.test(db=db)
pclassSRnat.high=pairedclass.test(db=db,var="SRnat")
pclassSRali.high=pairedclass.test(db=db,var="SRali")

### frequencies of significantly negative effects and thresholds
effects.classSR= summary.pc(pc=pclassSR, group="ALIEN");   mtext(text="Total richness", side=3, outer=T, line=3)
effects.classSRnat= summary.pc(pc=pclassSRnat, group="ALIEN") ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.classSRali= summary.pc(pc=pclassSRali, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=3)

effects.classSR.grass= summary.pc(pc=pclassSR.grass, group="ALIEN") ; mtext(text="Total richness", side=3, outer=T, line=3)
effects.classSRnat.grass= summary.pc(pc=pclassSRnat.grass, group="ALIEN") ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.classSRali.grass= summary.pc(pc=pclassSRali.grass, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=3)

effects.classSR.wood= summary.pc(pc=pclassSR.wood, group="ALIEN") ; mtext(text="total SR", side=3, outer=T, line=-2)
effects.classSRnat.wood= summary.pc(pc=pclassSRnat.wood, group="ALIEN") ; mtext(text="Native richness", side=3, outer=T, line=-2)
effects.classSRali.wood= summary.pc(pc=pclassSRali.wood, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=-2)

effects.classSR.low= summary.pc(pc=pclassSR.low, group="ALIEN") ; mtext(text="Total richness", side=3, outer=T, line=3)
effects.classSRnat.low= summary.pc(pc=pclassSRnat.low, group="ALIEN") ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.classSRali.low= summary.pc(pc=pclassSRali.low, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=3)


effects.classSR.high= summary.pc(pc=pclassSR.high, group="ALIEN") ; mtext(text="Total richness", side=3, outer=T, line=3)
effects.classSRnat.high= summary.pc(pc=pclassSRnat.high, group="ALIEN") ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.classSRali.high= summary.pc(pc=pclassSRali.high, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=3)


# boxplots for significant species
plot.pc(pclassSR, var= "SR")
plot.pc(pclassSRnat, var= "SRnat")
plot.pc(pclassSRnat, var= "SRali")
plot.pc(pclassSRnat.grass, var= "SRali", db = databp[databp$PlotName %in% grasslands,])

## proportion of thresholds per class
par (mfrow=c(1,3), mar=c(5,5,4,2))
m="natural"
threshold.spline(effects.classSR, m=m)
threshold.spline(effects.classSRnat, tit = "Native richness", m=m)
threshold.spline(effects.classSRali, tit = "Alien richness", m=m)

# density distrib
par (mfrow=c(1,3), mar=c(5,5,4,2))
m="natural"
threshold.spline(effects.classSR, index= "perc.thr",  m=m)
threshold.spline(effects.classSRnat, index= "perc.thr", tit = "Native richness", m=m)
threshold.spline(effects.classSRali, index= "perc.thr", tit = "Alien richness", m=m)

### thresholds summary
thresholds=data.frame(SR=pclassSR$threshold[rownames(pclassSR$threshold),1], SRnat=pclassSRnat$threshold[rownames(pclassSR$threshold),1],SRali=pclassSRali$threshold[rownames(pclassSR$threshold),1],
                      SR.grass=pclassSR.grass$threshold[rownames(pclassSR$threshold),1], SRnat.grass=pclassSRnat.grass$threshold[rownames(pclassSR$threshold),1],SRali.grass=pclassSRali.grass$threshold[rownames(pclassSR$threshold),1],
                      SR.wood=pclassSR.wood$threshold[rownames(pclassSR$threshold),1], SRnat.wood=pclassSRnat.wood$threshold[rownames(pclassSR$threshold),1],SRali.wood=pclassSRali.wood$threshold[rownames(pclassSR$threshold),1],
                      row.names=rownames(pclassSR$threshold))
thresholds$ALIEN= species[rownames(thresholds),"ALIEN"]
thresholds$sig=as.numeric(rowSums(thresholds[, c("SRnat", "SRali", "SR")], na.rm=T)>0)

tmp=thresholds
tmp[is.na(tmp)]=10
plot(jitter(tmp$SRali), jitter(tmp$SRnat), type='p', col=c("forestgreen","firebrick")[tmp$ALIEN + 1])
text(jitter(tmp$SRali), jitter(tmp$SRnat), label= rownames(tmp), offset=0, cex=0.5, col=c("forestgreen","firebrick")[tmp$ALIEN + 1])

library(FactoMineR)
pca1=PCA(thresholds[thresholds$sig==1,1:4], quali.sup=4)
x11()
par(mfrow=c(1,2))
plot(pca1, choix="ind", habillage=4) 
plot(pca1, choix="var", ann=F) 

# NMDS
library(vegan)
tmp=thresholds[thresholds$sig==1,1:3]
tmp = -(tmp-6)
tmp[is.na(tmp)]= 0
fit <- metaMDS (tmp, distance="gower", k=2,plot=T)
plot(fit, type="t")

k=kmeans(fit$points,2)

plot(fit$points, pch = c(21,10)[thresholds$ALIEN + 1], col = k$cluster)
arrows(rep(0,3),rep(0,3), fit$species[,1], fit$species[,2], col = "blue", pch=21)
text( x=fit$species[,1], y=fit$species[,2],label=rownames(fit$species),pos=1,  col = "blue")




######## ranking species according to their total impact
par(mfrow=c(1,3))
pc=pclassSR
plot(pc$threshold$nplots.impact , -pc$threshold$mean.sig.diff,
     main= "Total richness", col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 

pc=pclassSRnat
plot(pc$threshold$nplots.impact , -pc$threshold$mean.sig.diff,
     main= "Native richness", col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] )    
     
pc=pclassSRali
     plot(pc$threshold$nplots.impact , -pc$threshold$mean.sig.diff,
          main= "Alien richness", col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 

### Percent impact per sp vs. prevalence
par(mfrow=c(1,3))
pc=pclassSR
plot(pc$threshold$prevalence,pc$threshold$pplots.impact , ann=F,pch=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Total richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(pc$threshold$pplots.impact>0.5 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],pc$threshold$pplots.impact[tmp], label= rownames(pc$threshold)[tmp], adj=1 )
     
pc=pclassSRnat
plot(pc$threshold$prevalence,pc$threshold$pplots.impact , ann=F,pch=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Native richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(pc$threshold$pplots.impact>0.5 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],pc$threshold$pplots.impact[tmp], label= rownames(pc$threshold)[tmp], adj=1 )

pc=pclassSRali
plot(pc$threshold$prevalence,pc$threshold$pplots.impact , ann=F,pch=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Alien richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(pc$threshold$pplots.impact>0.5 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],pc$threshold$pplots.impact[tmp], label= rownames(pc$threshold)[tmp], adj=1 )


# grasslands
par(mfrow=c(1,3))
pc=pclassSR.grass
plot(pc$threshold$nplots.impact , -pc$threshold$mean.sig.diff,
     main= "Total richness", col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 

pc=pclassSRnat.grass
plot(pc$threshold$nplots.impact , -pc$threshold$mean.sig.diff,
     main= "Native richness", col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 

pc=pclassSRali.grass
plot(pc$threshold$nplots.impact , -pc$threshold$mean.sig.diff,
     main= "Alien richness", col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 

# woodlands
par(mfrow=c(1,3))
pc=pclassSR.wood
plot(pc$threshold$nplots.impact , -pc$threshold$mean.sig.diff,
     main= "Total richness", col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 

pc=pclassSRnat.wood
plot(pc$threshold$nplots.impact , -pc$threshold$mean.sig.diff,
     main= "Native richness", col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 

pc=pclassSRali.wood
plot(pc$threshold$nplots.impact , -pc$threshold$mean.sig.diff,
     main= "Alien richness", col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 



