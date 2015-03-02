# Figure Article thresholds of impact
## only results in real grasslands = first rank species is herbaceous + landcover is recorded as "grassland"

db <- databp[databp$PlotName %in% realgrasslands,]

### Preliminary figures :
### plot effects summary = proportion of significant negative effects and proportion of thresholds
plot.effect.summary (glmSR.sum$class.summary)
plot.effect.summary (glmSRnat.sum$class.summary)
plot.effect.summary (glmSRali.sum$class.summary)

#########  Details per significant species ##############
# raw observations for significant species
x11()
plot.glm(M=glmSR, var= "SR", db= db, sel.criteria = "th.exist", bst =T, ES = F, boxplots =F)
mtext(2,text = "total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat, var= "SRnat",  db= db,sel.criteria = "th.exist", bst =T, ES = F, boxplots =F)
mtext(2,text = "Native richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRali, var= "SRali",  db= db,sel.criteria = "th.exist", bst =T, ES = F, boxplots =F)
mtext(2,text = "Alien richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)


# boxplots for significant species
x11()
plot.glm(M=glmSR, var= "SR", db= db, sel.criteria = "th.exist", bst =T, ES = F)
mtext(2,text = "total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat, var= "SRnat",  db= db,sel.criteria = "th.exist", bst =T, ES = F)
mtext(2,text = "Native richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRali, var= "SRali",  db= db,sel.criteria = "th.exist", bst =T, ES = F)
mtext(2,text = "Alien richness", outer=T, line= 1, ES = F)
mtext(1,text = "abundance class", outer=T, line= 1)



# Effect size for significant species
x11()
plot.glm(M=glmSR, var= "SR", db= db, sel.criteria = "th.exist", bst =T, ES = T)
mtext(2,text = "effect size on total richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)

x11()
plot.glm(M = glmSRnat, var= "SRnat",  db= db,sel.criteria = "th.exist")
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


########## Impact spread #####

# Number of impacted plots vs. prevalence of species
plot.impact(x="prevalence", y="nb.plot.impact", square =T)
plot.impact(x="prevalence", y="prop.plot.impact", square =F)

plot.impact(x="prevalence", y=c("nb.plot.impact","/","prevalence"), square =F)

########### Impact size  #########
### Impact wtd mean size per sp vs. nb.impact
plot.impact(x="prop.plot.impact", y="wtd.mean.dif", square =F)
plot.impact(x="prop.plot.impact", y="prop.wt.mean.dif", square =F)

### Impact max size per sp vs. nb.impact
plot.impact(x="prop.plot.impact", y="mean.dif", square =F)
plot.impact(x="prop.plot.impact", y="prop.mean.dif", square =F)


### MAX Impact size per sp vs.nb.impact
plot.impact(x="prop.plot.impact", y="max.dif", square =F)
plot.impact(x="prop.plot.impact", y="prop.max.dif", square =F)

### threshold per sp vs. nb.impact
plot.impact(x="nb.plot.impact", y="th", square =F)

### threshold Impact size per sp vs. threshold
plot.impact(x="th", y="max.diff", square =F)
plot.impact(x="th", y="wtd.mean.diff", square =F)


### Ranking species by impact size * spread of impact
rank.impact(x="prop.wt.mean.dif", y="nb.plot.impact")
rank.impact(x="prop.wt.mean.dif", y="prop.plot.impact")
rank.impact(x="max.dif", y="prop.plot.impact")

rank.impact(x="index", y=NA)

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
