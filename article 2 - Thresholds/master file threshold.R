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
R <- 999
db <- databp[databp$PlotName %in% realgrasslands,]
min.occur <- 5
# 
# # perform analyses
#  system.time(glmSR <- glm.test(db = db,bootstrap = T, min.occur= min.occur, R=R, CI=0.95, drastic =F))
#  glmSRnat <- glm.test(db = db,var="SRnat",bootstrap = T,min.occur= min.occur,  R=R, CI=0.95, drastic =F)
#  glmSRali <- glm.test(db = db,var="SRali",bootstrap = T,min.occur= min.occur,  R=R, CI=0.95, drastic =F)
# 
# save(glmSR,glmSRnat,glmSRali, file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")

load(file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")

# ## GLMs with dominance index as covariate
# system.time(glmSR.dom <- glm.test(db = db,bootstrap = T, R=R,min.occur= min.occur,  covar ="dominance"))
# glmSRnat.dom <- glm.test(db = db,var="SRnat",bootstrap = T, R=R, min.occur= min.occur, covar ="dominance")
# glmSRali.dom <- glm.test(db = db,var="SRali",bootstrap = T, R=R,min.occur= min.occur,  covar ="dominance")
# save(glmSR.dom,glmSRnat.dom,glmSRali.dom, file = "saved Rdata/article 2 - threshold/booststrapped+dominance.glms.Rdata")


######## calculate proportional impacts (mean % of species lost)
glmSR$boot.thresh = add.prop(N = glmSR, var ="SR", data=db)
glmSRnat$boot.thresh = add.prop(N = glmSRnat, var ="SRnat", data=db)
glmSRali$boot.thresh = add.prop(N = glmSRali, var ="SRali", data=db)


### SUMMARY of frequencies of significantly negative effects and thresholds
glmSR.sum <-  summary.glmtest(M=glmSR, group="ALIEN", type="boot")
glmSRnat.sum <- summary.glmtest(M=glmSRnat, group="ALIEN", type="boot") 
glmSRali.sum <- summary.glmtest(M=glmSRali, group="ALIEN", type="boot") 

# save results
save.image("saved Rdata/article 2 - threshold/article threshold 1.0.Rdata")

### Statistical analyses

# create overall contingency tables


glm.boot.table <- lapply(list(Total.richness = glmSR, Native.richness = glmSRnat, Alien.richness =glmSRali), FUN= function (M) {
tbl <- as.table(t(rbind(none = as.numeric(table(na=is.na(M$boot.thresh$th), alien = rownames(M$boot.thresh)%in%aliens)[2,]),
                          table(M$boot.thresh$th, alien = rownames(M$boot.thresh)%in%aliens))))
return(tbl)
})

glm.table <- lapply(list(Total.richness = glmSR, Native.richness = glmSRnat, Alien.richness =glmSRali), FUN= function (M) {
  tbl <- as.table(t(rbind(none = as.numeric(table(na=is.na(M$thresh$th), alien = rownames(M$thresh)%in%aliens)[2,]),
                          table(M$thresh$th, alien = rownames(M$thresh)%in%aliens))))
  return(tbl)
})

## compare natives and aliens total threshold/impact occurrence :
tbl <- glmSRali.sum$overall.summary
tbl$nothr <- tbl$nb.sp - tbl$freq.thr 
tbl= t(tbl[,3:4])

chisq.test(tbl) ### significant for SRnat P = 0.02144

#GLM with binomial
# (thresh, no thresh) ~ class.abun + alien.status
tbl <- glmSRali.sum $class.summary

tbl$nothr <- tbl$n.target - tbl$freq.thr

fm1 <- cbind(freq.thr, nothr) ~ group
alienstatus_SR <- glm(fm1, data = tbl, family = binomial())
summary(alienstatus_SR)

fm1 <- cbind(freq.thr, nothr) ~ as.factor(class)
alienstatus_SR <- glm(fm1, data = tbl, family = binomial())
summary(alienstatus_SR)

fm2 <- cbind(freq.thr, nothr) ~ group + as.factor(class)
alienstatus_SR <- glm(fm2, data = tbl, family = binomial())
summary(alienstatus_SR)

## compare natives and aliens threshold occurrence per class:

#SR
tbl <- glmSR.sum $class.summary
chisq.test(unstack(tbl, form = nb.sp ~ group)) ### no difference in distrib of occurrence per class
chisq.test(unstack(tbl, form = n.target ~ group)) ### no difference in distrib of occurrence per class
t.test(n.target/nb.sp ~ group, data = tbl, paired = T) ### there is a difference in the proportion of target sp per class

chisq.test(unstack(tbl, form = freq.thr ~ group)) ### problem : too many zeros in native group
chisq.test(unstack(tbl, form = freq.negative ~ group)) ### diff in proportion of negative coeff
chisq.test(unstack(tbl, form = freq.positive ~ group)[1:4,]) ### no diff in proportion of positive coeff

t.test(freq.impact/nb.sp ~ group, data=tbl, paired = T)
t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T) 


#SRnat
tbl <- glmSRnat.sum $class.summary
chisq.test(unstack(tbl, form = n.target ~ group)) ### no difference in distrib of occurrence per class
t.test(n.target/nb.sp ~ group, data = tbl, paired = T) ### there is a difference in the proportion of target sp per class

chisq.test(unstack(tbl, form = freq.thr ~ group)[c(1,3,4),]) ### ns
chisq.test(unstack(tbl, form = freq.negative ~ group)) ### diff in proportion of negative coeff
chisq.test(unstack(tbl, form = freq.positive ~ group)[1:4,]) ### no diff in proportion of positive coeff

t.test(freq.impact/nb.sp ~ group, data=tbl, paired = T)  ### signif P= 0.01681, df=4, t=-3.9508
t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T)   ### signif P = 0.03048, df=4, t=-32808


#SRali
tbl <- glmSRnat.sum $class.summary
chisq.test(unstack(tbl, form = n.target ~ group)) ### no difference in distrib of occurrence per class
t.test(n.target/nb.sp ~ group, data = tbl, paired = T) ### there is a difference in the proportion of target sp per class

chisq.test(unstack(tbl, form = freq.thr ~ group)[c(1,3,4),]) ### ns
chisq.test(unstack(tbl, form = freq.negative ~ group)) ### diff in proportion of negative coeff
chisq.test(unstack(tbl, form = freq.positive ~ group)[1:4,]) ### no diff in proportion of positive coeff

t.test(freq.impact/nb.sp ~ group, data=tbl, paired = T)  ### signif P= 0.01681, df=4, t=-3.9508
t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T)   ### signif P = 0.03048, df=4, t=-32808


## other method : Independence tests : NOTHING SIGNIFICANT
tbl <- as.table(t(rbind(none = as.numeric(table(na=is.na(glmSRnat$boot.thresh$th), alien = rownames(glmSRnat$boot.thresh)%in%aliens)[2,]),
table(glmSRnat$boot.thresh$th, alien = rownames(glmSRnat$boot.thresh)%in%aliens))
 ))
library(coin)
cmh_test(tbl)
chisq_test(as.table( tbl))
