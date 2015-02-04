## Threshold master file

# setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R")

library(doBy)
library(vegan)
library(coin)

### import data
source('script/data/import BP species and environment data.R')

### taxonomy solving : 
# extract names, match them to DB, create a reference list of changes, update initial DB


# load functions
source('script/functions/p2star.R')
source('script/article 2 - Thresholds/glm test.R')


### threshold analysis using GLMs

## applying functions

# glmSR=glm.test()
# glmSRnat=glm.test(var="SRnat")
# glmSRali=glm.test(var="SRali")

glmSR.grass=glm.test(db=databp[databp$PlotName %in% realgrasslands,])
glmSRnat.grass=glm.test(db=databp[databp$PlotName %in% realgrasslands,],var="SRnat")
glmSRali.grass=glm.test(db=databp[databp$PlotName %in% realgrasslands,],var="SRali")

### frequencies of significantly negative effects and thresholds

# effects.glmSR= summary.glmtest(M=glmSR, group="ALIEN")
# effects.glmSRnat= summary.glmtest(M=glmSRnat, group="ALIEN") 
# effects.glmSRali= summary.glmtest(M=glmSRali, group="ALIEN") 

effects.glmSR.grass= summary.glmtest(M=glmSR.grass, group="ALIEN")
effects.glmSRnat.grass= summary.glmtest(M=glmSRnat.grass, group="ALIEN") 
effects.glmSRali.grass= summary.glmtest(M=glmSRali.grass, group="ALIEN") 

# save results
save.image("saved Rdata/article 2 - threshold/article threshold 1.0.Rdata")
