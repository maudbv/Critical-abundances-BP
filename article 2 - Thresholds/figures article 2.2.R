# Figure Article thresholds of impact
## only results in real grasslands = first rank species is herbaceous + landcover is recorded as "grassland"

db <- databp[databp$PlotName %in% realgrasslands,]

### Preliminary figures :
### plot effects summary = proportion of significant negative effects and proportion of thresholds
plot.effect.summary (effects.glmSR.grass)
plot.effect.summary (effects.glmSRnat.grass)
plot.effect.summary (effects.glmSRali.grass)

# boxplots for significant species
plot.glm(M=glmSR.grass, var= "SR", db= db, sel.criteria = "th.exist")
plot.glm(M = glmSRnat.grass, var= "SRnat",  db= db,sel.criteria = "th.exist")
plot.glm(M = glmSRali.grass, var= "SRali",  db= db,sel.criteria = "th.exist")

#### Frequency of impacts per class
par(mfcol=c(2,3), oma=c(3,6,4,2), mar=c(2,3,2,1))
thresh.freq(effects=effects.glmSR.grass, ylim=c(0,200))
thresh.freq(effects=effects.glmSRnat.grass, y=F, ylim=c(0,200))
thresh.freq(effects=effects.glmSRali.grass, y=F, ylim=c(0,200), leg=T)  

mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=0, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
mtext(text="Abundance class", side=1, outer=T, line=1)


########## Impact spread #####

# Number of impacted plots vs. prevalence of species
plot.impact(x="prevalence", y="nb.plot.impact", square =T)
plot.impact(x="prevalence", y="prop.plot.impact", square =F)


########### Impact size  #########
### Impact wtd mean size per sp vs. nb.impact
plot.impact(x="nb.plot.impact", y="wtd.mean.diff", square =F)

### Impact max size per sp vs. nb.impact
plot.impact(x="nb.plot.impact", y="max.diff", square =F)

### MAX Impact size per sp vs.nb.impact
plot.impact(x="nb.plot.impact", y="max.diff", square =F)

### threshold per sp vs. nb.impact
plot.impact(x="nb.plot.impact", y="th", square =F)

### threshold Impact size per sp vs. threshold
plot.impact(x="th", y="max.diff", square =F)
plot.impact(x="th", y="wtd.mean.diff", square =F)


### Ranking species by impact size * spread of impact
rank.impact(x="wtd.mean.diff", y="nb.plot.impact")
rank.impact(x="wtd.mean.diff", y="prop.plot.impact")
rank.impact(x="max.diff", y="nb.plot.impact")


##### impact index vs. time since introduciton
M <-glmSR.grass
y= (M$thresh$wtd.mean.diff * M$thresh$nb.plot.impact)
x=traitdata[rownames(M$thresh), "Nat_start"]
plot(x,y)
summary(lm(y~x)) # not signif

y=species$Sp.occurence
x=traitdata[rownames(species), "Nat_start"] #prevalence signif related to time since intro
plot(x,y)
summary(lm(y~x))
