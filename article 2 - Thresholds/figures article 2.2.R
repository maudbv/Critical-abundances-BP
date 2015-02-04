# Figure Article thresholds of impact



# boxplots for significant species
plot.glm(M=glmSR, var= "SR", db=databp, sel.criteria = "spear")
plot.glm(glmSRnat, var= "SRnat", sel.criteria = "spear")
plot.glm(glmSRali, var= "SRali", sel.criteria = "spear")
plot.glm(glmSR.grass, var= "SR", db = databp[databp$PlotName %in% grasslands,], sel.criteria = "spear")



# Fig 1

### all plots
par(mfcol=c(2,3), oma=c(3,7,4,2), mar=c(2,3,2,1))
thresh.freq( ylim=c(0,300))
thresh.freq(effects=effects.glmSRnat, y=F, ylim=c(0,300))
thresh.freq(effects=effects.glmSRali, y=F, ylim=c(0,300), leg=T)

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=5, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=0, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)



##  GRasslands

par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.freq(effects=effects.glmSR.grass, ylim=c(0,200))
thresh.freq(effects=effects.glmSRnat.grass, y=F, ylim=c(0,200))
thresh.freq(effects=effects.glmSRali.grass, y=F, ylim=c(0,200), leg=T)  

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=0, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)

# representing proportion of threshold
par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.prop(effects=effects.glmSR.grass)
thresh.prop(effects=effects.glmSRnat.grass, y=F)
thresh.prop(effects=effects.glmSRali.grass, y=F)  

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=-1, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)

### compare to LSD
glmsp = which(rowSums(glmSR.grass$P<=0.05, na.rm=T)>0)
length(glmsp)
sum(!is.na(glmSR.grass$glms$th))

lsdsp = which(abs(rowSums(LSDclassSR.grass$leastdiff, na.rm=T))>0)
length(lsdsp)
sum(!is.na(LSDclassSR.grass$threshold$threshold))

which(!is.na(LSDclassSR.grass$threshold$threshold))  %in%  which(!is.na(glmSR.grass$glms$th))

plot(unlist(glmSR.grass$z), unlist(LSDclassSR.grass$leastdiff + LSDclassSR.grass$diff), xlab="GLM z-score", ylab = "LSD - diff")
abline(v=0,h=0)


row.names(grass.effect.summary[grass.effect.summary$total.effect>0,]) %in% names(glmsp)
row.names(grass.effect.summary[grass.effect.summary$total.effect>0,]) %in% rownames(glmSR.grass$glms)[which(!is.na(glmSR.grass$glms$th))]



##  GRasslands

par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.freq(effects=effects.LSDclassSR.grass, ylim=c(0,200))
thresh.freq(effects=effects.LSDclassSRnat.grass, y=F, ylim=c(0,200))
thresh.freq(effects=effects.LSDclassSRali.grass, y=F, ylim=c(0,200), leg=T)  

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=-1, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)

# representing proportion of threshold
par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.prop(effects=effects.LSDclassSR.grass)
thresh.prop(effects=effects.LSDclassSRnat.grass, y=F)
thresh.prop(effects=effects.LSDclassSRali.grass, y=F)  


# table of all significant species in grasslands
cols=c("threshold","pplots.impact","max.sig.diff","thr.diff")          
grassland.thr=data.frame(SR=cbind(LSDclassSR.grass$threshold[, c("prevalence", cols)], round(LSDclassSR.grass$k$diff, 4)),
                         SRnat=cbind(LSDclassSRnat.grass$threshold[,  cols], round(LSDclassSRnat.grass$k$diff, 4)),
                         SRali=cbind(LSDclassSRali.grass$threshold[,  cols], round(LSDclassSRali.grass$k$diff, 4)))
grassland.thr = grassland.thr[which(rowSums(grassland.thr[, c("SR.pplots.impact", "SRnat.pplots.impact","SRali.pplots.impact")]>0, na.rm=T)>0),]
grassland.thr$ALIEN=species[rownames(grassland.thr), "ALIEN"]

write.csv(as.matrix(grassland.thr), file="table thresholds grassland.csv")

##  Woodlands
par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.freq(effects=effects.LSDclassSR.wood, ylim=c(0,200))
thresh.freq(effects=effects.LSDclassSRnat.wood, y=F, ylim=c(0,200))
thresh.freq(effects=effects.LSDclassSRali.wood, y=F, ylim=c(0,200), leg=T)  

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=-1, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)


# table of all significant species in woodlands

woodland.thr=data.frame(SR=cbind(LSDclassSR.wood$threshold, LSDclassSR.wood$k),
                         SRnat=cbind(LSDclassSRnat.wood$threshold, LSDclassSRnat.wood$k),
                         SRali=cbind(LSDclassSRali.wood$threshold, LSDclassSRali.wood$k))
woodland.thr = woodland.thr[which(rowSums(woodland.thr[, c("SR.nplots.impact", "SRnat.nplots.impact","SRali.nplots.impact")]>0, na.rm=T)>0),]
woodland.thr$ALIEN=species[rownames(woodland.thr), "ALIEN"]

write.csv(as.matrix(woodland.thr), file="table thresholds woodland.csv")

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
