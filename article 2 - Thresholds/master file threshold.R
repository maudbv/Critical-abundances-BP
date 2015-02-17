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

# set bootstrap sample size
R=999
db=databp[databp$PlotName %in% realgrasslands,]

# perform analyses
#  system.time(glmSR.grass <- glm.test(db = db,bootstrap = T, R=R))
#  glmSRnat.grass <- glm.test(db = db,var="SRnat",bootstrap = T, R=R)
#  glmSRali.grass <- glm.test(db = db,var="SRali",bootstrap = T, R=R)
# save(glmSR.grass,glmSRnat.grass,glmSRali.grass, file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")
# 
 load(file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")
# 
# system.time(glmSR.grass <- glm.test(db = db,bootstrap = F, R=R))
# glmSRnat.grass <- glm.test(db = db,var="SRnat",bootstrap = F, R=R)
# glmSRali.grass <- glm.test(db = db,var="SRali",bootstrap = F, R=R)
# save(glmSR.grass,glmSRnat.grass,glmSRali.grass, file = "saved Rdata/article 2 - threshold/non-booststrapped.glms.Rdata")

# system.time(glmSR.grass <- glm.test(db = db,bootstrap = T, R=R, covar ="dominance"))
# glmSRnat.grass <- glm.test(db = db,var="SRnat",bootstrap = T, R=R, covar ="dominance")
# glmSRali.grass <- glm.test(db = db,var="SRali",bootstrap = T, R=R, covar ="dominance")
# save(glmSR.grass,glmSRnat.grass,glmSRali.grass, file = "saved Rdata/article 2 - threshold/booststrapped+dominance.glms.Rdata")

### frequencies of significantly negative effects and thresholds
effects.glmSR.grass <-  summary.glmtest(M=glmSR.grass, group="ALIEN", type="boot")
effects.glmSRnat.grass <- summary.glmtest(M=glmSRnat.grass, group="ALIEN", type="boot") 
effects.glmSRali.grass <- summary.glmtest(M=glmSRali.grass, group="ALIEN", type="boot") 

# save results
save.image("saved Rdata/article 2 - threshold/article threshold 1.0.Rdata")
load("saved Rdata/article 2 - threshold/article threshold 1.0.Rdata")
