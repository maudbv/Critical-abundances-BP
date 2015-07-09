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
                            glmSRali.sum.th$class.summary )))

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
         ES = T, threshold = threshold, sp = sp.target, panel = c(3,5))
mtext(2,text = "effect size on total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat.overall, var= "SRnat",  db= db,sel.criteria = "sp.target", type= "overall.boot",
         ES = T, threshold = threshold, sp = sp.target, panel = c(4,3),ylim = c(-3,3))
mtext(2,text = "effect size on Native richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRali.overall, var= "SRali",  db= db,sel.criteria = "sp.target", type= "overall.boot",
         ES = T, threshold = threshold, sp = sp.target, ylim = c(-0.5,0.5), panel = c(5,3))
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

### Graphs of delta gamma vs delta alpha 
#community change plots
## species order
cm <- comm[(rownames(comm) %in% realgrasslands),]
cm <- cm[,colnames(cm)%in%natives]
cm <- ceiling(cm>0)
cm <- cm[,colSums(cm)>0]
ord <- colnames(cm)[order(colSums(cm))]
ord <- ord[ord %in% natives]

x11()
par( mar=c(1,3,1,0))

gams <-  rowSums(divpart.nat.perm$obs[,c("shared","lost", "gained")])
alphas <- as.numeric(impact.SRnat[rownames(divpart.nat.perm$obs),]$SRo)

dg<- (divpart.nat.perm$obs$above.gamm - divpart.nat.perm$obs$below.gamm)/gams
da<- (divpart.nat.perm$obs$above.alpha - divpart.nat.perm$obs$below.alpha) / alphas

gams <-  rowSums(divpart.ali.perm$obs[,c("shared","lost", "gained")])
alphas <- as.numeric(impact.SRali[rownames(divpart.nat.perm$obs),]$SRo)

dga<- (divpart.ali.perm$obs$above.gamm - divpart.ali.perm$obs$below.gamm)/gams
daa<- (divpart.ali.perm$obs$above.alpha - divpart.ali.perm$obs$below.alpha) / alphas

dg<- (divpart.nat.perm$P$above.gamm -0.5)*2
da<- (divpart.nat.perm$P$above.alpha - 0.5)*2
dga<- (divpart.ali.perm$P$above.gamm -0.5)*2
daa<- (divpart.ali.perm$P$above.alpha - 0.5)*2


cgam <- c("darkgrey", "black")[ (divpart.nat.perm$P$above.gamm<=0.05 | divpart.nat.perm$P$above.gamm>=0.95) + 1]

plot(da, dg, xlim =c(-1,1),ylim =c(-1,1) , pch = 20, col=cgam )
points(daa, dga, xlim =c(-1,1),ylim =c(-1,1) , pch = "+", col=cgam )
segments (da, dg, daa, dga,col=cgam)
text(x = daa + c(0.07, 0.07, 0.15,+0.15,0.15, 0.07, 0.07,0.07, 0.07, 0.07, 0.07),
     y = dga +c(0.05, 0.05, 0.05,+0.1,0.05, 0.05, 0.05,0.05, 0.05, 0.05, 0.05),
     label = rownames(divpart.nat.perm$P) , col=cgam, cex =0.7 )

abline(h=0, v=0, col="grey")
