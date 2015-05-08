# Figure Article thresholds of impact
## only results in real grasslands = first rank species is herbaceous + landcover is recorded as "grassland"

db <- databp[databp$PlotName %in% realgrasslands,]
sp.target <-  rownames(glmSR.overall$impact.spread) [
  which(  rownames(glmSRnat.overall$impact.spread) %in% aliens &
         ( !is.na(glmSRnat.overall$impact.spread$th.CI) |
        !is.na(glmSRali.overall$impact.spread$th.CI)))]

save(sp.target, file = "saved Rdata/article 2 - threshold/sp.target.Rdata")

### Preliminary figures :
### plot effects summary = proportion of significant negative effects and proportion of thresholds
plot.effect.summary ()
plot.effect.summary (effects =
                       list(glmSR.sum.th$class.summary,
                            glmSRnat.sum.th$class.summary,
                            glmSRali.sum.th$class.summary ))

plot.effect.summary.freq ()
#########  Details per significant species ##############
threshold = "th.CI"

# raw observations for significant species
x11()
plot.glm(M=glmSR.overall, var= "SR", db= db, type= "overall.boot", ES = F,  panels = c(4,6),boxplots =F, threshold = threshold)
mtext(2,text = "total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat.overall, var= "SRnat",  db= db,sel.criteria = "th.exist", type= "overall.boot", ES = F, boxplots =F, threshold = threshold)
mtext(2,text = "total richness", outer=T, line= 1)
mtext(2,text = "Native richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRali.overall, var= "SRali",  db= db,sel.criteria = "th.exist", type= "overall.boot", ES = F, boxplots =F, threshold = threshold)
mtext(2,text = "total richness", outer=T, line= 1)
mtext(2,text = "Alien richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)


# boxplots for significant species
x11()
plot.glm(M=glmSR.overall, var= "SR", db= db, sel.criteria = "sp.target", type= "overall.boot", 
         ES = F, threshold = threshold, sp = sp.target)
mtext(2,text = "total richness", outer=T, line= 1)
mtext(2,text = "total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat.overall, var= "SRnat",  db= db,sel.criteria = "sp.target", type= "overall.boot",
         ES = F, threshold = threshold, sp = sp.target)
mtext(2,text = "Native richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRali.overall, var= "SRali",  db= db,sel.criteria = "sp.target", type= "overall.boot",
         ES = F, threshold = threshold, sp = sp.target)
mtext(2,text = "Alien richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)



# Effect size for significant species
x11()
plot.glm(M=glmSR.overall, var= "SR", db= db,sel.criteria = "sp.target", type= "overall.boot",
         ES = T, threshold = threshold, sp = sp.target)
mtext(2,text = "effect size on total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat.overall, var= "SRnat",  db= db,sel.criteria = "sp.target", type= "overall.boot",
         ES = T, threshold = threshold, sp = sp.target)
mtext(2,text = "effect size on Native richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRali.overall, var= "SRnat",  db= db,sel.criteria = "sp.target", type= "overall.boot",
         ES = T, threshold = threshold, sp = sp.target)
mtext(2,text = "effect size on Alien richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)




# 
# ### frequency with bootstrap variance
# x11()
# 
# tab <- apply(glmSR.overall$crit.valsCI[,], 2, FUN=function(x) {
#   x <- factor(x,levels= c("2","3","4","5","6"))
#   f = table(x)
# })
# 
# sd = apply(tab,1,  FUN= sd, na.rm=T)
# barplot(tab[match(2:6, rownames(tab)),1], ylim= c(0,max(tab) +2))
# segments(1:4, tab[,1] - sd, 1:4, tab[,1] + sd)
# mtext(2,text = "Frequency critical value", outer=F, line= 3)
# mtext(1,text = "abundance class", outer=F, line= 3)

########## Impact spread #####

# Number of impacted plots vs. prevalence of species

par(mfrow=c(1,1),las = 1, mar=c(2,2,2,2) ,oma=c(2,2,0,0), cex = 0.9)
for (i in 1:1) {
    M <- list( glmSRnat.overall,glmSRali.overall) [[i]]
       sp <- which(rownames(M$impact.spread)%in%aliens & M$impact.spread$prop.plot.impact >0)  ## selects only aliens
       M=lapply(M, FUN=function(x) x=x[sp,])

    M$impact.spread$prop.plot.impact = M$impact.spread$prop.plot.impact*100
    
plot(prop.plot.impact ~ prevalence , data = M$impact.spread, xlim=c(-20, 800),
     ylim=c(0, 85),ann=F,pch=20, bg = "black", cex.lab = 0.8)

dpt <- as.matrix(dist(M$impact.spread[,c("prop.plot.impact","prevalence")]))<30
dpt[upper.tri(dpt)] <- NA
diag(dpt)<- NA
offsets <- rep(0, length(sp))
names(offsets) <- sp
offsets[which(rowSums(dpt, na.rm = T) >= 1)] <- 1
offsets[which(colSums(dpt, na.rm = T)>=1)] <- -1


text(M$impact.spread$prevalence + offsets*60, M$impact.spread$prop.plot.impact +abs(offsets)*2 ,
     label = rownames(M$impact.spread), cex = 0.7, pos = 3, offset = 0.3)

segments(M$impact.spread$prevalence , M$impact.spread$prop.plot.impact,
         M$impact.spread$prevalence + offsets*60, M$impact.spread$prop.plot.impact +abs(offsets)*2)

# mtext(3, text = c("Native richness", "Alien richness")[i])

}
mtext(1, text ="Number of occurrences",outer=T, line=0.5)
mtext(2, text = "% > critical abundance", outer = T,las = 0, line =1)



########### Impact size  #########
### Impact wtd mean size per sp vs. nb.impact

plot.impact(x="prop.plot.impact", y="prop.wtd.mean.dif", square =F,xlab="prop. plots > critical abundance", ylab="average decrease in %SR")

### MAX Impact size per sp vs.nb.impact
plot.impact(x="prop.plot.impact", y="prop.max.dif", square =F,xlab="prop. plots > critical abundance", ylab="max decrease in %SR")

### threshold per sp vs. nb.impact
plot.impact(x="nb.plot.impact", y="th", square =F)

### representing impact index
x="index"

i=2
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
    
    levelplot(z2~ S * P, col.regions= colorRampPalette(c("beige" , "firebrick")),
              main = "impact index", xlim= c(0,1), ylim=c(0,1))
    trellis.focus("panel", 1, 1, highlight=FALSE)
    lpoints(M$thresh$prop.wtd.mean.dif,M$thresh$prop.plot.impact  , pch=3,col="black", cex=1)
    
    abs <- M$thresh$prop.wtd.mean.dif +0.01
    ord <- M$thresh$prop.plot.impact + runif(length( M$thresh$prop.plot.impact), 0, +0.01)
    ltext(abs, ord, label= as.character(rownames(M$thresh)), 
          col="black", cex=0.8, adj =0 )
    trellis.unfocus()
    



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
