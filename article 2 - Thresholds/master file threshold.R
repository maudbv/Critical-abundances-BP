## Threshold master file

setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R")

library(doBy)
library(vegan)

### import data
source('script/data/import BP species and environment data.R')

### taxonomy solving : 
# extract names, match them to DB, create a reference list of changes, update initial DB

# save.image(file="saved Rdata/Banks peninsula dataset.Rdata")


### plot subsampling for thresholds : only real realgrasslands

### import functions
source('script/article 2 - Thresholds/glm test.R')


### threshold analysis using GLMs

## applying functions

glmSR=glm.test()
glmSRnat=glm.test(var="SRnat")
glmSRali=glm.test(var="SRali")

glmSR.grass=glm.test(db=databp[databp$PlotName %in% realgrasslands,])
glmSRnat.grass=glm.test(db=databp[databp$PlotName %in% realgrasslands,],var="SRnat")
glmSRali.grass=glm.test(db=databp[databp$PlotName %in% realgrasslands,],var="SRali")

### frequencies of significantly negative effects and thresholds
effects.glmSR= summary.glmtest(M=glmSR, group="ALIEN", graph=F);   mtext(text="Total richness", side=3, outer=T, line=3)
effects.glmSRnat= summary.glmtest(M=glmSRnat, group="ALIEN", graph=F) ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.glmSRali= summary.glmtest(M=glmSRali, group="ALIEN", graph=F) ; mtext(text="Alien richness", side=3, outer=T, line=3)

effects.glmSR.grass= summary.glmtest(M=glmSR.grass, group="ALIEN", graph=F) ; mtext(text="Total richness", side=3, outer=T, line=3)
effects.glmSRnat.grass= summary.glmtest(M=glmSRnat.grass, group="ALIEN", graph=F) ;  mtext(text="Native richness", side=3, outer=T, line=3)
effects.glmSRali.grass= summary.glmtest(M=glmSRali.grass, group="ALIEN", graph=F) ; mtext(text="Alien richness", side=3, outer=T, line=3)
