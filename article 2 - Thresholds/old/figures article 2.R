# Figure Article thresholds of impact


# Fig 1

thresh.prop=function(effects=effects.classSR, data=species, y=T,ylim=c(0,0.35)) {
    #graphical representation
    barplot(effects$prop.thr[effects$group=="ALIEN:0"], names.arg= paste("c",2:6, sep=""), 
            col="grey", ylim=ylim, las=1 )
    if (y==T) mtext(side=2, text="Prop. threshold effects",line=3.4, cex=0.8)
    barplot( effects$prop.thr[effects$group=="ALIEN:1"], names.arg= paste("c",2:6, sep=""),
             col="black" , ylim=ylim , las=1)
    if (y==T) mtext(side=2, text="Prop. threshold effects", line=3.4, cex=0.8)
    
    print(wilcox.test(effects$prop.thr[effects$group=="ALIEN:1"],effects$prop.thr[effects$group=="ALIEN:0"], paired=T))
    print(friedman.test(cbind(effects$prop.thr[effects$group=="ALIEN:1"],effects$prop.thr[effects$group=="ALIEN:0"])))
        
}

# thresh.freq=function(effects=effects.classSR, data=species, y=T,ylim=c(0,100), leg=F) {
#   #graphical representation
#   tmp=effects
#   tmp$total=tmp$nb.sp - tmp$freq.impact
#   tmp$nonthr= tmp$freq.impact - tmp$freq.thr
#   
#   x=t(as.matrix(tmp[tmp$group=="ALIEN:0",c("freq.thr", "nonthr", "total")]))
#   barplot(x, names.arg= paste("c",2:6, sep=""),ylim=ylim, las=1, 
#           col=c("black", "grey", "white"))
#   if (y==T) mtext(side=2, text="Number of species",line=3.4, cex=0.8)
#  
#   x=t(as.matrix(tmp[tmp$group=="ALIEN:1",c("freq.thr", "nonthr", "total")]))
#   barplot(x, names.arg= paste("c",2:6, sep=""), ylim=ylim, las=1, 
#           col=c("black", "grey", "white"))
#   if (y==T) mtext(side=2, text="Number of species", line=3.4, cex=0.8)
#   if (leg) {
#     legend("topright", bty="n", legend=c("threshold", "impact", "no impact"), fill=c("black", "grey", "white"), cex=0.9)
#   }
# chisq.test(t(cbind(effects$freq.thr[effects$group=="ALIEN:1"],effects$freq.thr[effects$group=="ALIEN:0"])))      
# }
# 
# 
# par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
# thresh.prop()
# thresh.prop(effects=effects.classSRnat, y=F)
# thresh.prop(effects=effects.classSRali, y=F)  
# 
# mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
# mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=-1, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
# mtext(text="Abundance class", side=1, outer=T, line=1)


### logarithmic scale

thresh.freq=function(effects=effects.classSR, data=species, y=T,ylim=c(0,100), leg=F) {
  tmp=effects
  # Native targets
  x=t(as.matrix(tmp[tmp$group=="ALIEN:0",c("freq.thr", "freq.impact", "nb.sp")]))
#    x[x!=0] = log(x[x!=0])
  x = log(x)
  barplot(x[3,]+1, names= paste("c",2:6, sep=""),las=1, ylim=c(0, max(log(ylim + 1))),
          col=c( "white"), yaxt="n")
  barplot(x[2,]+1, col=c( "grey"), names="",axes=F,add=T)
  barplot(x[1,]+1, col=c("black"),  names="",axes=F, add=T)
  axis(side=2, at=log(c(1,2,5,10,25,50,100, 200))+ 1, label=c(1,2,5,10,25,50,100,200), las=1)
  box(bty="l")
  if (y==T) mtext(side=2, text="Number of species", line=3.4, cex=0.8)

  if (leg) {
    legend("topright", bty="n", legend=c("threshold", "impact", "no impact"), fill=c("black", "grey", "white"), cex=0.9)
  }
  
  ### Alien targets
  x=t(as.matrix(tmp[tmp$group=="ALIEN:1",c("freq.thr", "freq.impact", "nb.sp")]))
  #    x[x!=0] = log(x[x!=0])
  x = log(x)
  barplot(x[3,]+1, names= paste("c",2:6, sep=""),las=1, ylim=c(0, max(log(ylim + 1))),
          col=c( "white"), yaxt="n")
  barplot(x[2,]+1, col=c( "grey"), names="",axes=F,add=T)
  barplot(x[1,]+1, col=c("black"),  names="",axes=F, add=T)
  axis(side=2, at=log(c(1,2,5,10,25,50,100, 200))+ 1, label=c(1,2,5,10,25,50,100,200), las=1)
  box(bty="l")
  if (y==T) mtext(side=2, text="Number of species", line=3.4, cex=0.8)
  
    chisq.test(t(cbind(effects$freq.thr[effects$group=="ALIEN:1"],effects$freq.thr[effects$group=="ALIEN:0"])))      
}


par(mfcol=c(2,3), oma=c(3,7,4,2), mar=c(2,3,2,1))
thresh.freq( ylim=c(0,300))
thresh.freq(effects=effects.classSRnat, y=F, ylim=c(0,300))
thresh.freq(effects=effects.classSRali, y=F, ylim=c(0,300), leg=T)

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=5, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=0, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)



##  GRasslands

par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.freq(effects=effects.classSR.grass, ylim=c(0,100))
thresh.freq(effects=effects.classSRnat.grass, y=F, ylim=c(0,100))
thresh.freq(effects=effects.classSRali.grass, y=F, ylim=c(0,100), leg=T)  

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=-1, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)


# table of all significant species in grasslands
cols=c("threshold","pplots.impact","max.sig.diff","thr.diff")          
grassland.thr=data.frame(SR=cbind(pclassSR.grass$threshold[, c("prevalence", cols)], round(pclassSR.grass$k$p.val, 4)),
                         SRnat=cbind(pclassSRnat.grass$threshold[,  cols], round(pclassSRnat.grass$k$p.val, 4)),
                         SRali=cbind(pclassSRali.grass$threshold[,  cols], round(pclassSRali.grass$k$p.val, 4)))
grassland.thr = grassland.thr[which(rowSums(grassland.thr[, c("SR.pplots.impact", "SRnat.pplots.impact","SRali.pplots.impact")]>0, na.rm=T)>0),]
grassland.thr$ALIEN=species[rownames(grassland.thr), "ALIEN"]

write.csv(as.matrix(grassland.thr), file="table thresholds grassland.csv")

##  Woodlands
par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.freq(effects=effects.classSR.wood, ylim=c(0,200))
thresh.freq(effects=effects.classSRnat.wood, y=F, ylim=c(0,200))
thresh.freq(effects=effects.classSRali.wood, y=F, ylim=c(0,200), leg=T)  

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=-1, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)


# table of all significant species in woodlands

woodland.thr=data.frame(SR=cbind(pclassSR.wood$threshold, pclassSR.wood$k),
                         SRnat=cbind(pclassSRnat.wood$threshold, pclassSRnat.wood$k),
                         SRali=cbind(pclassSRali.wood$threshold, pclassSRali.wood$k))
woodland.thr = woodland.thr[which(rowSums(woodland.thr[, c("SR.nplots.impact", "SRnat.nplots.impact","SRali.nplots.impact")]>0, na.rm=T)>0),]
woodland.thr$ALIEN=species[rownames(woodland.thr), "ALIEN"]

write.csv(as.matrix(woodland.thr), file="table thresholds woodland.csv")

# Fig 2 
### Percent impact per sp vs. prevalence

par(mfrow=c(1,3), cex=0.8)

pc=pclassSR 
sp=which(rownames(pclassSR$threshold)%in%aliens)

pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$prevalence,pc$threshold$pplots.impact , ann=F,pch=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Total richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(pc$threshold$pplots.impact>0.5 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],pc$threshold$pplots.impact[tmp], 
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7 )

pc=pclassSRnat
sp=which(rownames(pclassSR$threshold)%in%aliens)
pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$prevalence,pc$threshold$pplots.impact , ann=F,pch=20, ylim=c(0,1), cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Native richness",ylab="% potential impact", xlab="Species Prevalence")
abline(h=0.5, lty="dashed")
tmp=which(pc$threshold$pplots.impact>0.5 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],pc$threshold$pplots.impact[tmp], 
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7 )

pc=pclassSRali
sp=which(rownames(pclassSR$threshold)%in%aliens)
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
sp=which(rownames(pclassSR$threshold)%in%aliens)
  
par(mfrow=c(1,3), cex=0.8)
  
pc=pclassSR 
pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$nplots.impact , -pc$threshold$wtd.mean.sig.diff, ann=F,pch=20,cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Total richness",ylab="impact size", xlab="frequency of potential impact")
abline(h=3, lty="dashed")
tmp=which(-pc$threshold$wtd.mean.sig.diff>3 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],-pc$threshold$wtd.mean.sig.diff[tmp],
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)

pc=pclassSRnat
pc=lapply(pc, FUN=function(x) x=x[sp,])
plot(pc$threshold$nplots.impact , -pc$threshold$wtd.mean.sig.diff, ann=F,pch=20, cex=1.2,
     col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
title( main= "Native richness",ylab="impact size", xlab="frequency of potential impact")
abline(h=3, lty="dashed")
tmp=which(-pc$threshold$wtd.mean.sig.diff>3 & pc$threshold$prevalence>600)
text(pc$threshold$prevalence[tmp],-pc$threshold$wtd.mean.sig.diff[tmp],
     label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)

pc=pclassSRali
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
  sp=which(rownames(pclassSR$threshold)%in%aliens)
  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=pclassSR 
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$nplots.impact , -pc$threshold$max.sig.diff, ann=F,pch=20,cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Total richness",ylab="Max impact size", xlab="frequency of potential impact")
  abline(h=3, lty="dashed")
  tmp=which(-pc$threshold$max.sig.diff>3 & pc$threshold$prevalence>600)
  text(pc$threshold$prevalence[tmp],-pc$threshold$max.sig.diff[tmp],
       label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)
  
  pc=pclassSRnat
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$nplots.impact , -pc$threshold$max.sig.diff, ann=F,pch=20, cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Native richness",ylab="Max impact size", xlab="frequency of potential impact")
  abline(h=3, lty="dashed")
  tmp=which(-pc$threshold$max.sig.diff>3 & pc$threshold$prevalence>600)
  text(pc$threshold$prevalence[tmp],-pc$threshold$max.sig.diff[tmp],
       label= rownames(pc$threshold)[tmp], adj=1, cex=0.7)
  
  pc=pclassSRali
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
  sp=which(rownames(pclassSR$threshold)%in%aliens)
  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=pclassSR 
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$wtd.mean.sig.diff, ann=F,pch=20,cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Total richness",ylab="impact size", xlab="threshold")
    
  pc=pclassSRnat
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$wtd.mean.sig.diff, ann=F,pch=20, cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Native richness",ylab="impact size", xlab="threshold")
    
  pc=pclassSRali
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$wtd.mean.sig.diff,  ann=F,pch=20,  cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Alien richness",ylab="impact size", xlab="threshold")
  
}) ()


### MAX Impact size per sp vs. threshold

(function() {
  
  # select only alien targets
  sp=which(rownames(pclassSR$threshold)%in%aliens)
  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=pclassSR 
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$max.sig.diff, ann=F,pch=20,cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Total richness",ylab="Max impact size", xlab="threshold")
  
  pc=pclassSRnat
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$max.sig.diff, ann=F,pch=20, cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Native richness",ylab="Max impact size", xlab="threshold")
  
  pc=pclassSRali
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$max.sig.diff,  ann=F,pch=20,  cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Alien richness",ylab="Max impact size", xlab="threshold")
  
}) ()


### threshold Impact size per sp vs. threshold

(function() {
  
  # select only alien targets
 
  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=pclassSR 
  sp=which(rownames(pclassSR$threshold)%in%aliens)
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$thr.diff, ann=F,pch=20,cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Total richness" , ylab="Threshold effect", xlab="threshold abundance")
  
  pc=pclassSRnat
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$thr.diff, ann=F,pch=20, cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Native richness",ylab="Threshold effect", xlab="threshold abundance")
  
  pc=pclassSRali
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  plot(pc$threshold$threshold , -pc$threshold$thr.diff,  ann=F,pch=20,  cex=1.2,
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1] ) 
  title( main= "Alien richness", ylab="Threshold effect", xlab="threshold abundance")
  
}) ()



### Ranking species per % of plots over threshold

(function() {
  
  # select only alien targets

  
  par(mfrow=c(1,3), cex=0.8)
  
  pc=pclassSR 
  sp=which(rownames(pclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact, decreasing=T),]
  plot(1:length(sp), pc$threshold$pplots.impact ,ann=F,type="h",cex=1.2, ylim=c(0,1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Total richness",ylab="% potential impact")
  
  pc=pclassSRnat
  sp=which(rownames(pclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact, decreasing=T),]
  plot(1:length(sp), pc$threshold$pplots.impact ,ann=F,type="h",cex=1.2, ylim=c(0,1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Native richness",ylab="% potential impact")
  
  pc=pclassSRali  
  sp=which(rownames(pclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
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
  
  pc=pclassSR 
  sp=which(rownames(pclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) , decreasing=T),]
  x=pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) 
  plot(1:length(sp),x ,ann=F,type="h",cex=1.2,
       ylim= c(0,max(x)+1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Total richness",ylab="Potential Impact index")
  
  pc=pclassSRnat
  sp=which(rownames(pclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
  pc=lapply(pc, FUN=function(x) x=x[sp,])
  pc$threshold=pc$threshold[order(pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) , decreasing=T),]
  x=pc$threshold$pplots.impact * (-pc$threshold$wtd.mean.sig.diff) 
  plot(1:length(sp),x ,ann=F,type="h",cex=1.2,
       ylim= c(0,max(x)+1),
       col=c("forestgreen", "firebrick")[ rownames(pc$threshold)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(pc$threshold), las=3, cex.axis=0.7)
  title( main= "Native richness",ylab="Potential Impact index")
  
  pc=pclassSRali  
  sp=which(rownames(pclassSR$threshold)%in%aliens & ! is.na(pc$threshold$pplots.impact))
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
y=- (pclassSR$threshold$thr.diff * pclassSR$threshold$nplots.impact)
y=pclassSR$threshold$nplots.impact
y=- pclassSR$threshold$thr.diff
x=traitdata[rownames(pclassSR$threshold), "Nat_start"]

plot(x,y)
summary(lm(y~x)) # not signif

y=species$Sp.occurence
x=traitdata[rownames(species), "Nat_start"] #prevalence signif related to time since intro
plot(x,y)
summary(lm(y~x))
