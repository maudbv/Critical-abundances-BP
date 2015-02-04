# Figure Article thresholds of impact
## only results in real grasslands = first rank species is herbaceous + landcover is recorded as "grassland"

# boxplots for significant species
plot.glm(M=glmSR.grass, var= "SR", db=databp, sel.criteria = "dev1")
plot.glm(glmSRnat.grass, var= "SRnat", sel.criteria = "dev1")
plot.glm(glmSRali.grass, var= "SRali", sel.criteria = "dev1")

# Fig 1
par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.freq(effects=effects.glmSR.grass, ylim=c(0,200))
thresh.freq(effects=effects.glmSRnat.grass, y=F, ylim=c(0,200))
thresh.freq(effects=effects.glmSRali.grass, y=F, ylim=c(0,200), leg=T)  

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=0, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)

### plot effects summary = proportion of significant negative effects and proportion of thresholds
plot.effect.summary (effects.glmSR.grass)
plot.effect.summary (effects.glmSRnat.grass)
plot.effect.summary (effects.glmSRali.grass)

# Fig 2 
### Percent impact per sp vs. prevalence

par(mfrow=c(1,3), cex=0.8)

pc=LSDclassSR 
sp=which(rownames(LSDclassSR$threshold)%in%aliens)

pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$prevalence,pc$threshold$pplots.impact , ann=F,pch=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Total richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(pc$threshold$pplots.impact>0.5 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],pc$threshold$pplots.impact[tmp], 
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7 )

pc=LSDclassSRnat
sp=which(rownames(LSDclassSR$threshold)%in%aliens)
pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$prevalence,pc$threshold$pplots.impact , ann=F,pch=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Native richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(pc$threshold$pplots.impact>0.5 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],pc$threshold$pplots.impact[tmp], 
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7 )

pc=LSDclassSRali
sp=which(rownames(LSDclassSR$threshold)%in%aliens)
pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$prevalence,pc$threshold$pplots.impact , ann=F,pch=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Alien richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(pc$threshold$pplots.impact>0.5 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],pc$threshold$pplots.impact[tmp], 
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7 )

### Impact size per sp vs. prevalence

(function() {

  # select only alien targets
sp=which(rownames(LSDclassSR$threshold)%in%aliens)
  
par(mfrow=c(1,3), cex=0.8)
  
pc=LSDclassSR 
pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$nplots.impact , -pc$threshold$wtd.mean.sig.diff, ann=F,pch=20,cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Total richness",ylab="impact size", xlab="frequency of potential impact")
abline(h=3, lty="dashed")
tmp=which(-pc$threshold$wtd.mean.sig.diff>3 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],-pc$threshold$wtd.mean.sig.diff[tmp],
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)

pc=LSDclassSRnat
pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$nplots.impact , -pc$threshold$wtd.mean.sig.diff, ann=F,pch=20, cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Native richness",ylab="impact size", xlab="frequency of potential impact")
abline(h=3, lty="dashed")
tmp=which(-pc$threshold$wtd.mean.sig.diff>3 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],-pc$threshold$wtd.mean.sig.diff[tmp],
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)

pc=LSDclassSRali
  pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$nplots.impact , -pc$threshold$wtd.mean.sig.diff,  ann=F,pch=20,  cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Alien richness",ylab="impact size", xlab="frequency of potential impact")
abline(h=3, lty="dashed")
tmp=which(-pc$threshold$wtd.mean.sig.diff>3 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],-pc$threshold$wtd.mean.sig.diff[tmp],
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)

}) ()

### MAX Impact size per sp vs. prevalence

(function() {
  
  # select only alien targets
  sp=which(rownames(LSDclassSR$threshold)%in%aliens)
  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=LSDclassSR 
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$nplots.impact , -pc$threshold$max.sig.diff, ann=F,pch=20,cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Total richness",ylab="Max impact size", xlab="frequency of potential impact")
  abline(h=3, lty="dashed")
  tmp=which(-pc$threshold$max.sig.diff>3 & pc$threshold$prevalence>600)
  text(pc$threshold$prevalence[tmp],-pc$threshold$max.sig.diff[tmp],
       label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)
  
  pc=LSDclassSRnat
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$nplots.impact , -pc$threshold$max.sig.diff, ann=F,pch=20, cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Native richness",ylab="Max impact size", xlab="frequency of potential impact")
  abline(h=3, lty="dashed")
  tmp=which(-pc$threshold$max.sig.diff>3 & pc$threshold$prevalence>600)
  text(pc$threshold$prevalence[tmp],-pc$threshold$max.sig.diff[tmp],
       label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)
  
  pc=LSDclassSRali
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$nplots.impact , -pc$threshold$max.sig.diff,  ann=F,pch=20,  cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Alien richness",ylab="Max impact size", xlab="frequency of potential impact")
  abline(h=3, lty="dashed")
  tmp=which(-pc$threshold$max.sig.diff>3 & pc$threshold$prevalence>600)
  text(pc$threshold$prevalence[tmp],-pc$threshold$max.sig.diff[tmp],
       label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)
  
}) ()

### Impact size per sp vs. threshold

(function() {
  
  # select only alien targets
  sp=which(rownames(LSDclassSR$threshold)%in%aliens)
  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=LSDclassSR 
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$wtd.mean.sig.diff, ann=F,pch=20,cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Total richness",ylab="impact size", xlab="threshold")
    
  pc=LSDclassSRnat
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$wtd.mean.sig.diff, ann=F,pch=20, cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Native richness",ylab="impact size", xlab="threshold")
    
  pc=LSDclassSRali
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$wtd.mean.sig.diff,  ann=F,pch=20,  cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Alien richness",ylab="impact size", xlab="threshold")
  
}) ()


### MAX Impact size per sp vs. threshold

(function() {
  
  # select only alien targets
  sp=which(rownames(LSDclassSR$threshold)%in%aliens)
  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=LSDclassSR 
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$max.sig.diff, ann=F,pch=20,cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Total richness",ylab="Max impact size", xlab="threshold")
  
  pc=LSDclassSRnat
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$max.sig.diff, ann=F,pch=20, cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Native richness",ylab="Max impact size", xlab="threshold")
  
  pc=LSDclassSRali
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$max.sig.diff,  ann=F,pch=20,  cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Alien richness",ylab="Max impact size", xlab="threshold")
  
}) ()


### threshold Impact size per sp vs. threshold

(function() {
  
  # select only alien targets
 
  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=LSDclassSR 
  sp=which(rownames(LSDclassSR$threshold)%in%aliens)
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$thr.diff, ann=F,pch=20,cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Total richness" , ylab="Threshold effect", xlab="threshold abundance")
  
  pc=LSDclassSRnat
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$thr.diff, ann=F,pch=20, cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Native richness",ylab="Threshold effect", xlab="threshold abundance")
  
  pc=LSDclassSRali
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$thr.diff,  ann=F,pch=20,  cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Alien richness", ylab="Threshold effect", xlab="threshold abundance")
  
}) ()



### Ranking species per % of plots over threshold

(function() {
  
  # select only alien targets

  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=LSDclassSR 
  sp=which(rownames(LSDclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact, decreasing=T),]
  plot(1:length(sp), pc$threshold$pplots.impact ,ann=F,type="h",cex=1.2, ylim=c(0,1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Total richness",ylab="% potential impact")
  
  pc=LSDclassSRnat
  sp=which(rownames(LSDclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact, decreasing=T),]
  plot(1:length(sp), pc$threshold$pplots.impact ,ann=F,type="h",cex=1.2, ylim=c(0,1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Native richness",ylab="% potential impact")
  
  pc=LSDclassSRali  
  sp=which(rownames(LSDclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact, decreasing=T),]
  plot(1:length(sp), pc$threshold$pplots.impact ,ann=F,type="h",cex=1.2, ylim=c(0,1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Alien richness",ylab="% potential impact")
  
}) ()

### Ranking species by impact size * frequenc of impact
(function() {

  par(mfrow=c(1,3), cex=0.8)
  
  pc=LSDclassSR 
  sp=which(rownames(LSDclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) , decreasing=T),]
  x=pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) 
  plot(1:length(sp),x ,ann=F,type="h",cex=1.2,
       ylim= c(0,max(x)+1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Total richness",ylab="Potential Impact index")
  
  pc=LSDclassSRnat
  sp=which(rownames(LSDclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) , decreasing=T),]
  x=pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) 
  plot(1:length(sp),x ,ann=F,type="h",cex=1.2,
       ylim= c(0,max(x)+1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Native richness",ylab="Potential Impact index")
  
  pc=LSDclassSRali  
  sp=which(rownames(LSDclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) , decreasing=T),]
  x=pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) 
  plot(1:length(sp),x ,ann=F,type="h",cex=1.2,
       ylim= c(0,max(x)+1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Alien richness",ylab="Potential Impact index")

}) ()




##### impact index vs. time since introduciton
y=- (LSDclassSR$threshold$thr.diff * LSDclassSR$threshold$nplots.impact)
y=LSDclassSR$threshold$nplots.impact
y=- LSDclassSR$threshold$thr.diff
x=traitdata[rownames(LSDclassSR$threshold), "Nat_start"]

plot(x,y)
summary(lm(y~x)) # not signif

y=species$Sp.occurence
x=traitdata[rownames(species), "Nat_start"] #prevalence signif related to time since intro
plot(x,y)
summary(lm(y~x))
