# Pairwise abundance classes tests

# source('script/data/import fede data.R')
# source('script/functions/function paired-class SR.R')
# Calculate t-tests between pairs of abundance classes :
LSDclassSR=LSDclass.test()
LSDclassSRnat=LSDclass.test(var="SRnat")
LSDclassSRali=LSDclass.test(var="SRali")

LSDclassSR.grass=LSDclass.test(db=databp[databp$PlotName %in% grasslands,])
LSDclassSRnat.grass=LSDclass.test(db=databp[databp$PlotName %in% grasslands,],var="SRnat")
LSDclassSRali.grass=LSDclass.test(db=databp[databp$PlotName %in% grasslands,],var="SRali")

LSDclassSR.wood=LSDclass.test(db=databp[databp$PlotName %in% woodlands,])
LSDclassSRnat.wood=LSDclass.test(db=databp[databp$PlotName %in% woodlands,],var="SRnat")
LSDclassSRali.wood=LSDclass.test(db=databp[databp$PlotName %in% woodlands,],var="SRali")

## divided by low and high grasslands

db=databp[databp$PlotName %in% grasslands & databp$DEM_10 < 400,]
length(unique(db$PlotName))

LSDclassSR.low=LSDclass.test(db=db)
LSDclassSRnat.low=LSDclass.test(db=db,var="SRnat")
LSDclassSRali.low=LSDclass.test(db=db,var="SRali")

db=databp[databp$PlotName %in% grasslands & databp$DEM_10 >= 400,]
length(unique(db$PlotName))
LSDclassSR.high=LSDclass.test(db=db)
LSDclassSRnat.high=LSDclass.test(db=db,var="SRnat")
LSDclassSRali.high=LSDclass.test(db=db,var="SRali")

### frequencies of significantly negative effects and thresholds
effects.LSDclassSR= summary.lsd(M=LSDclassSR, group="ALIEN");   mtext(text="Total richness", side=3, outer=T, line=3)
effects.LSDclassSRnat= summary.lsd(M=LSDclassSRnat, group="ALIEN") ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.LSDclassSRali= summary.lsd(M=LSDclassSRali, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=3)

effects.LSDclassSR.grass= summary.lsd(M=LSDclassSR.grass, group="ALIEN") ; mtext(text="Total richness", side=3, outer=T, line=3)
effects.LSDclassSRnat.grass= summary.lsd(M=LSDclassSRnat.grass, group="ALIEN") ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.LSDclassSRali.grass= summary.lsd(M=LSDclassSRali.grass, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=3)

effects.LSDclassSR.wood= summary.lsd(M=LSDclassSR.wood, group="ALIEN") ; mtext(text="total SR", side=3, outer=T, line=-2)
effects.LSDclassSRnat.wood= summary.lsd(M=LSDclassSRnat.wood, group="ALIEN") ; mtext(text="Native richness", side=3, outer=T, line=-2)
effects.LSDclassSRali.wood= summary.lsd(M=LSDclassSRali.wood, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=-2)

effects.LSDclassSR.low= summary.lsd(M=LSDclassSR.low, group="ALIEN") ; mtext(text="Total richness", side=3, outer=T, line=3)
effects.LSDclassSRnat.low= summary.lsd(M=LSDclassSRnat.low, group="ALIEN") ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.LSDclassSRali.low= summary.lsd(M=LSDclassSRali.low, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=3)


effects.LSDclassSR.high= summary.lsd(M=LSDclassSR.high, group="ALIEN") ; mtext(text="Total richness", side=3, outer=T, line=3)
effects.LSDclassSRnat.high= summary.lsd(M=LSDclassSRnat.high, group="ALIEN") ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.LSDclassSRali.high= summary.lsd(M=LSDclassSRali.high, group="ALIEN") ; mtext(text="Alien richness", side=3, outer=T, line=3)


# boxplots for significant species
plot.lsd(LSDclassSR, var= "SR")
plot.lsd(LSDclassSRnat, var= "SRnat")
plot.lsd(LSDclassSRali, var= "SRali")
plot.lsd(LSDclassSR.grass, var= "SR", db = databp[databp$PlotName %in% grasslands,])

## proportion of thresholds per class
par (mfrow=c(1,3), mar=c(5,5,4,2))
m="natural"
threshold.spline(effects.LSDclassSR, m=m)
threshold.spline(effects.LSDclassSRnat, tit = "Native richness", m=m)
threshold.spline(effects.LSDclassSRali, tit = "Alien richness", m=m)




######## ranking species according to their total impact
par(mfrow=c(1,3))
lsd=LSDclassSR
plot(lsd$threshold$nplots.impact , -lsd$threshold$mean.sig.diff,
     main= "Total richness", col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 

lsd=LSDclassSRnat
plot(lsd$threshold$nplots.impact , -lsd$threshold$mean.sig.diff,
     main= "Native richness", col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] )    
     
lsd=LSDclassSRali
     plot(lsd$threshold$nplots.impact , -lsd$threshold$mean.sig.diff,
          main= "Alien richness", col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 

### Percent impact per sp vs. prevalence
par(mfrow=c(1,3))
lsd=LSDclassSR
plot(lsd$threshold$prevalence,lsd$threshold$pplots.impact , ann=F,lsdh=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 
title( main= "Total richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(lsd$threshold$pplots.impact>0.5 & lsd$threshold$prevalence>600)
text(lsd$threshold$prevalence[tmp],lsd$threshold$pplots.impact[tmp], label= rownames(lsd$threshold)[tmp], adj=1 )
     
lsd=LSDclassSRnat
plot(lsd$threshold$prevalence,lsd$threshold$pplots.impact , ann=F,lsdh=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 
title( main= "Native richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(lsd$threshold$pplots.impact>0.5 & lsd$threshold$prevalence>600)
text(lsd$threshold$prevalence[tmp],lsd$threshold$pplots.impact[tmp], label= rownames(lsd$threshold)[tmp], adj=1 )

lsd=LSDclassSRali
plot(lsd$threshold$prevalence,lsd$threshold$pplots.impact , ann=F,lsdh=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 
title( main= "Alien richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(lsd$threshold$pplots.impact>0.5 & lsd$threshold$prevalence>600)
text(lsd$threshold$prevalence[tmp],lsd$threshold$pplots.impact[tmp], label= rownames(lsd$threshold)[tmp], adj=1 )


# grasslands
par(mfrow=c(1,3))
lsd=LSDclassSR.grass
plot(lsd$threshold$nplots.impact , -lsd$threshold$mean.sig.diff,
     main= "Total richness", col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 

lsd=LSDclassSRnat.grass
plot(lsd$threshold$nplots.impact , -lsd$threshold$mean.sig.diff,
     main= "Native richness", col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 

lsd=LSDclassSRali.grass
plot(lsd$threshold$nplots.impact , -lsd$threshold$mean.sig.diff,
     main= "Alien richness", col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 

# woodlands
par(mfrow=c(1,3))
lsd=LSDclassSR.wood
plot(lsd$threshold$nplots.impact , -lsd$threshold$mean.sig.diff,
     main= "Total richness", col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 

lsd=LSDclassSRnat.wood
plot(lsd$threshold$nplots.impact , -lsd$threshold$mean.sig.diff,
     main= "Native richness", col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 

lsd=LSDclassSRali.wood
plot(lsd$threshold$nplots.impact , -lsd$threshold$mean.sig.diff,
     main= "Alien richness", col=c("forestgreen", "firebrick")[ rownames(lsd$threshold)%in% aliens +1] ) 



