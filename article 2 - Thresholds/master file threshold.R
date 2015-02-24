## Threshold master file

setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
# setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R")

library(doBy)
library(vegan)
library(coin)

# load("saved Rdata/article 2 - threshold/article threshold 1.0.Rdata")

### import data
# load("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/saved Rdata/article 3 - trait and phylo/save article 3.Rdata.RData")

source('script/data/import BP species and environment data.R', encoding = "native.enc")
source('script/data/import trait data.R', encoding = "native.enc")

### taxonomy solving : 
# extract names, match them to DB, create a reference list of changes, update initial DB


# load functions
source('script/functions/p2star.R')
source('script/article 2 - Thresholds/glm test.R')
source('script/article 2 - Thresholds/plotting functions.R')


### threshold analysis using GLMs

# set bootstrap sample size
# R <- 9999
# db <- databp[databp$PlotName %in% realgrasslands,]
# min.occur <- 3
# # perform analyses
#  system.time(glmSR <- glm.test(db = db,bootstrap = T, min.occur= min.occur, R=R))
#  glmSRnat <- glm.test(db = db,var="SRnat",bootstrap = T,min.occur= min.occur,  R=R)
#  glmSRali <- glm.test(db = db,var="SRali",bootstrap = T,min.occur= min.occur,  R=R)
# save(glmSR,glmSRnat,glmSRali, file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")

load(file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")

# ## GLMs with dominance index as covariate
# system.time(glmSR.dom <- glm.test(db = db,bootstrap = T, R=R,min.occur= min.occur,  covar ="dominance"))
# glmSRnat.dom <- glm.test(db = db,var="SRnat",bootstrap = T, R=R, min.occur= min.occur, covar ="dominance")
# glmSRali.dom <- glm.test(db = db,var="SRali",bootstrap = T, R=R,min.occur= min.occur,  covar ="dominance")
# save(glmSR.dom,glmSRnat.dom,glmSRali.dom, file = "saved Rdata/article 2 - threshold/booststrapped+dominance.glms.Rdata")

### frequencies of significantly negative effects and thresholds
glmSR.sum <-  summary.glmtest(M=glmSR, group="ALIEN", type="boot")
glmSRnat.sum <- summary.glmtest(M=glmSRnat, group="ALIEN", type="boot") 
glmSRali.sum <- summary.glmtest(M=glmSRali, group="ALIEN", type="boot") 

# save results
save.image("saved Rdata/article 2 - threshold/article threshold 1.0.Rdata")

### Statistical analyses

# export contingency tables
glm.table <- lapply(list(Total.richness = glmSR, Native.richness = glmSRnat, Alien.richness =glmSRali), FUN= function (M) {
tbl <- as.table(t(rbind(none = as.numeric(table(na=is.na(M$boot.thresh$th), alien = rownames(M$boot.thresh)%in%aliens)[2,]),
                          table(M$boot.thresh$th, alien = rownames(M$boot.thresh)%in%aliens))))
return(tbl)
})

## compare natives and aliens total threshold/impact occurrence :
tbl <- glmSRali.sum$overall.summary
tbl$nothr <- tbl$nb.sp - tbl$freq.thr 
tbl= t(tbl[,3:4])

chisq.test(tbl) 


## compare natives and aliens threshold occurrence per class:

tbl <- glmSR.sum $class.summary
chisq.test(unstack(tbl, form = nb.sp ~ group)) ### no difference in distrib of occurrence per class

# compare proportion of threshold in each class between alien and native targets
t.test(n.target ~ group, data = tbl, paired = T) 
t.test(freq.impact/nb.sp ~ group, data=tbl, paired = T)

t.test(freq.impact/n.target ~ group, data = tbl, paired = T) 
t.test(freq.thr/n.target ~ group, data = tbl, paired =T) 

## other method : Independence tests : nope
tbl <- as.table(t(rbind(none = as.numeric(table(na=is.na(glmSRnat$boot.thresh$th), alien = rownames(glmSRnat$boot.thresh)%in%aliens)[2,]),
table(glmSRnat$boot.thresh$th, alien = rownames(glmSRnat$boot.thresh)%in%aliens))
 ))

cmh_test(tbl)
chisq_test(as.table( tbl[,-1]))

