## Threshold master file

setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
# setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R")

library(doBy)
library(vegan)
library(coin)

### import data
# load("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/saved Rdata/article 3 - trait and phylo/save article 3.Rdata.RData")

source('script/data/import BP species and environment data.R', encoding = "native.enc")
source('script/data/import trait data.R', encoding = "native.enc")

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

# set bootstrap sample size
R=999
db=databp[databp$PlotName %in% realgrasslands,]

# perform analyses
system.time(glmSR.grass <- glm.test(db = db,bootstrap = T, R=R))

"Error during wrapup: number of items to replace is not a multiple of replacement length"

glmSRnat.grass <- glm.test(db = db,var="SRnat",bootstrap = T, R=R)
glmSRali.grass <- glm.test(db = db,var="SRali",bootstrap = T, R=R)

### frequencies of significantly negative effects and thresholds

# effects.glmSR <- summary.glmtest(M=glmSR, group="ALIEN")
# effects.glmSRnat= summary.glmtest(M=glmSRnat, group="ALIEN") 
# effects.glmSRali= summary.glmtest(M=glmSRali, group="ALIEN") 

effects.glmSR.grass <-  summary.glmtest(M=glmSR.grass, group="ALIEN")
effects.glmSRnat.grass <- summary.glmtest(M=glmSRnat.grass, group="ALIEN") 
effects.glmSRali.grass <- summary.glmtest(M=glmSRali.grass, group="ALIEN") 

# save results
save.image("saved Rdata/article 2 - threshold/article threshold 1.0.Rdata")
