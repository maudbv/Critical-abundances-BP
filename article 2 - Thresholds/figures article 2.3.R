# Figure Article thresholds of impact
## only results in real grasslands = first rank species is herbaceous + landcover is recorded as "grassland"

db <- databp[databp$PlotName %in% realgrasslands,]

### Preliminary figures :
### plot effects summary = proportion of significant negative effects and proportion of thresholds
plot.effect.summary ()

#########  Details per significant species ##############
# raw observations for significant species
x11()
plot.glm(M=glmSR.overall, var= "SR", db= db, type= "overall.boot", ES = F,  panels = c(5,4),boxplots =F)
mtext(2,text = "total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat, var= "SRnat",  db= db,sel.criteria = "th.exist", type= "overall.boot", ES = F, boxplots =F)
mtext(2,text = "Native richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRali, var= "SRali",  db= db,sel.criteria = "th.exist", type= "overall.boot", ES = F, boxplots =F)
mtext(2,text = "Alien richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)


# boxplots for significant species
x11()
plot.glm(M=glmSR, var= "SR", db= db, sel.criteria = "th.exist", type= "overall.boot", ES = F)
mtext(2,text = "total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat, var= "SRnat",  db= db,sel.criteria = "th.exist", type= "overall.boot", ES = F)
mtext(2,text = "Native richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRali, var= "SRali",  db= db,sel.criteria = "th.exist", type= "overall.boot", ES = F)
mtext(2,text = "Alien richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)



# Effect size for significant species
x11()
plot.glm(M=glmSR, var= "SR", db= db, sel.criteria = "th.exist", type= "overall.boot", panel= c(5,4),ES = T)
mtext(2,text = "effect size on total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat, var= "SRnat",  db= db,sel.criteria = "th.exist",  panels = c(5,4))
mtext(2,text = "effect size on Native richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRali, var= "SRali",  db= db,sel.criteria = "th.exist")
mtext(2,text = "effect size on Alien richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)



#### ### Frequency of impacts per class ###############
x11()

par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.freq(effects=glmSR.sum$class.summary, ylim=c(0,200))
thresh.freq(effects=glmSRnat.sum$class.summary, y=F, ylim=c(0,200))
thresh.freq(effects=glmSRali.sum$class.summary, y=F, ylim=c(0,200), leg=T)  

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=0, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)

### frequency with bootstrap variance
x11()
tab <- apply(glmSR.overall$crit.vals, 2, FUN=function(x) {
  x <- factor(x,levels= c("2","3","4","5","6"))
  f = table(x)
})

sd = apply(tab,1,  FUN= sd, na.rm=T)
barplot(tab[match(2:6, rownames(tab)),1], ylim= c(0,max(tab) +2))
segments(1:4, tab[,1] - sd, 1:4, tab[,1] + sd)

########## Impact spread #####

# Number of impacted plots vs. prevalence of species
plot.impact(x="prevalence", y="prop.plot.impact", square =F, xlab="total occurrences", ylab="proportion of plots > critical abundance")

########### Impact size  #########
### Impact wtd mean size per sp vs. nb.impact

plot.impact(x="prop.plot.impact", y="prop.wtd.mean.dif", square =F,xlab="prop. plots > critical abundance", ylab="average decrease in %SR")

### MAX Impact size per sp vs.nb.impact
plot.impact(x="prop.plot.impact", y="prop.max.dif", square =F,xlab="prop. plots > critical abundance", ylab="max decrease in %SR")

### threshold per sp vs. nb.impact
plot.impact(x="nb.plot.impact", y="th", square =F)

### representing impact index
x="index"

  for (i in 1:3) {
    M <- list(glmSR, glmSRnat, glmSRali) [[i]]
    M$thresh <- M$boot.thresh
    
    sp <- which(rownames(M$thresh)%in%aliens & ( M$thresh$nb.plot.impact>0))
    M=lapply(M, FUN=function(x) x=x[sp,])
    M$thresh$impact <-  M$thresh[,x]
    
    o  <- order(M$thresh$impact, decreasing=T)
    M=lapply(M, FUN=function(x) x=x[o,])
    
    t = seq(0,1, 0.001)
    S <- matrix(t, length(t),length(t), byrow=T)
    P <- matrix(t, length(t),length(t),byrow=F)
    z2 <- sqrt(S*P)
    
    library(lattice)
    
    x11() 
    
    levelplot(z2~ S * P, col.regions= colorRampPalette(c("beige" , "firebrick", "purple")),
              main = "impact index", xlim= c(0,1), ylim=c(0,1))
    trellis.focus("panel", 1, 1, highlight=FALSE)
    lpoints(M$thresh$prop.wtd.mean.dif,M$thresh$prop.plot.impact  , pch=3,col="black", cex=1)
    
    abs <- M$thresh$prop.wtd.mean.dif +0.01
    ord <- M$thresh$prop.plot.impact + runif(length( M$thresh$prop.plot.impact), 0, +0.01)
    ltext(abs, ord, label= as.character(rownames(M$thresh)), 
          col="black", cex=0.8, adj =0 )
    trellis.unfocus()
    
  }


### Ranking species by impact size * spread of impact
x="index"
x11()
par(mfrow=c(1,3), cex=0.8, oma=c(2,2,2,2), mar=c(4,2,1,1))
for (i in 1:3) {
  M <- list(glmSR, glmSRnat, glmSRali) [[i]]
  M$thresh <- M$boot.thresh
  
  sp <- which(rownames(M$thresh)%in%aliens & ( M$thresh$nb.plot.impact>0))
  M=lapply(M, FUN=function(x) x=x[sp,])
  M$thresh$impact <-  M$thresh[,x]
  
  o  <- order(M$thresh$impact, decreasing=T)
  M=lapply(M, FUN=function(x) x=x[o,])
  

  plot(1:length(sp), M$thresh$impact,ann=F,type="h",cex=1.2,
       ylim= c(0,0.5),
       col=c("forestgreen", "firebrick")[ rownames(M$thresh)%in% aliens +1], xaxt="n") 
  axis(side=1, at=1:length(sp), label=rownames(M$thresh), las=3, cex.axis=0.7)
  title( main= c("Total richness", "Native Richness","Alien Richness")[i])
}
mtext(2, text = paste("impact index") , outer= T, line=1)



##### impact index vs. time since introduciton
M <-glmSR
y= (M$boot.thresh$prop.wt.mean.diff * M$boot.thresh$nb.plot.impact)
x=traitdata[rownames(M$thresh), "Nat_start"]
plot(x,y)
summary(lm(y~x)) # not signif

y=species$Sp.occurence
x=traitdata[rownames(species), "Nat_start"] #prevalence signif related to time since intro
plot(x,y)
summary(lm(y~x))
