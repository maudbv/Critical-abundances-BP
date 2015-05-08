## Statistical analyses


# summarize thresholds


thresh.summary <- data.frame(SR = glmSR.overall$impact.spread[,c("th", "th.CI")], SRnat = glmSRnat.overall$impact.spread[,c("th", "th.CI")], SRali = glmSRali.overall$impact.spread[,c("th", "th.CI")], 
                             status = c( "NATIVE","ALIEN")[rownames(glmSR.overall$impact.spread) %in% aliens +1])

write.csv(thresh.summary, file = "threshold.summary.csv")

# create overall contingency tables
threshold = "th.CI"
glm.table <- lapply(list(Total.richness = glmSR.overall, Native.richness = glmSRnat.overall, Alien.richness =glmSRali.overall), FUN= function (M) {
  tbl <- as.table(t(rbind(none = as.numeric(table(na=is.na(M$impact.spread[, threshold]), alien = rownames(M$impact.spread)%in%aliens)[2,]),
                          table(M$impact.spread[, threshold], alien = rownames(M$impact.spread)%in%aliens))))
  return(tbl)
})



## compare natives and aliens total threshold/impact occurrence :
tbl<- glmSR.sum$overall.summary
tbl$nothr <- tbl$nb.sp - tbl$freq.thr 
tbl= t(tbl[,3:4])

chisq.test(tbl)

tbl <- glmSRnat.sum$overall.summary
tbl$nothr <- tbl$nb.sp - tbl$freq.thr 
tbl= t(tbl[,3:4])

chisq.test(tbl) 



tbl <- glmSRali.sum$overall.summary
tbl$nothr <- tbl$nb.sp - tbl$freq.thr 
tbl= t(tbl[,3:4])

chisq.test(tbl)

### counting number of species in common for three components of richness :
rownames(glmSRnat.overall$impact.spread)[which(!is.na(glmSRnat.overall$impact.spread$th) & !is.na(glmSR.overall$impact.spread$th))]
rownames(glmSRnat.overall$impact.spread)[which(is.na(glmSRnat.overall$impact.spread$th) & is.na(glmSR.overall$impact.spread$th))]
rownames(glmSRnat.overall$impact.spread)[which(!is.na(glmSRnat.overall$impact.spread$th) & !is.na(glmSRali.overall$impact.spread$th))]


# VGLM using VGAM package
library(VGAM)
# (thresh, no thresh) ~ class.abun + alien.status

tbl <- glmSR.sum $class.summary
vglm.SR<- vglm(freq.thr ~ class + group, family = poissonff,
                 data = tbl, trace = TRUE)
summary(vglm.SR)

tbl <- glmSRnat.sum $class.summary
vglm.SRnat<- vglm(freq.thr ~ class + group, family = poissonff,
               data = tbl, trace = TRUE)
summary(vglm.SRnat)

tbl <- glmSRali.sum $class.summary
vglm.SRali<- vglm(freq.thr ~ class + group, family = poissonff,
                  data = tbl, trace = TRUE)
summary(vglm.SRali)

#GLM with binomial
# (thresh, no thresh) ~ class.abun + alien.status
tbl <- glmSRnat.sum $class.summary
tbl$nothr <- tbl$nb.sp - tbl$freq.thr

fm1 <- cbind(freq.thr, nothr) ~ group
alienstatus_SR <- glm(fm1, data = tbl, family = binomial())
summary(alienstatus_SR)

fm1 <- cbind(freq.thr, nothr) ~ as.factor(class)
alienstatus_SR <- glm(fm1, data = tbl, family = binomial())
summary(alienstatus_SR)

fm2 <- cbind(freq.thr, nothr) ~ group + as.factor(class)
alienstatus_SR <- glm(fm2, data = tbl, family = binomial())
summary(alienstatus_SR)

## compare natives and aliens threshold occurrence per class: Mann whitney U two-sample test

tbl <- glmSR.sum $class.summary
wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$freq.th[tbl$group =="ALIEN:1"],
            alternative ="two.sided", paired=T)

tbl <- glmSRnat.sum $class.summary
wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$freq.th[tbl$group =="ALIEN:1"],
            alternative ="two.sided", paired=T)

tbl <- glmSRali.sum $class.summary
wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$freq.th[tbl$group =="ALIEN:1"],
            alternative ="two.sided", paired=T)

## compare natives and aliens threshold occurrence per class:mean rannks method
tbl <- glmSR.sum $class.summary
wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$freq.th[tbl$group =="ALIEN:1"],
            alternative ="two.sided", paired=T)

tbl <- glmSRnat.sum $class.summary
wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$freq.th[tbl$group =="ALIEN:1"],
            alternative ="two.sided", paired=T)

tbl <- glmSRali.sum $class.summary
wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$freq.th[tbl$group =="ALIEN:1"],
            alternative ="two.sided", paired=T)


## compare natives and aliens threshold occurrence per class:

#SR
tbl <- glmSR.sum $class.summary
t.test(freq.negative/nb.sp ~ group, data = tbl, paired = T) 
t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T) 

round(tbl$freq.thr/tbl$nb.sp,4)


#SRnat
tbl <- glmSRnat.sum $class.summary
t.test(freq.negative/nb.sp ~ group, data = tbl, paired = T) ### there is a difference in the proportion of target sp per class
t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T)   ### signif P = 0.03048, df=4, t=-32808
t.test(freq.thr/freq.negative ~ group, data = tbl, paired =T)   
round(tbl$freq.thr/tbl$nb.sp , 4)


wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"]/tbl$freq.negative[tbl$group =="ALIEN:0"],
            tbl$freq.th[tbl$group =="ALIEN:1"]/tbl$freq.negative[tbl$group =="ALIEN:1"],
            alternative ="two.sided", paired=T)


#SRali
tbl <- glmSRali.sum $class.summary
t.test(freq.negative/nb.sp ~ group, data = tbl, paired = T) 
t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T)   
t.test(freq.thr/freq.negative ~ group, data = tbl, paired =T)   



round(tbl$freq.thr/tbl$nb.sp, 4)

# 
# 
# chisq.test(unstack(tbl, form = n.target ~ group)) ### no difference in distrib of occurrence per class
# t.test(n.target/nb.sp ~ group, data = tbl, paired = T) ### there is a difference in the proportion of target sp per class
# t.test(freq.negative/nb.sp ~ group, data = tbl, paired = T) ### there is a difference in the proportion of target sp per class
# 
# chisq.test(unstack(tbl, form = freq.thr ~ group)[c(1,3,4),]) ### ns
# chisq.test(unstack(tbl, form = freq.negative ~ group)) ### diff in proportion of negative coeff
# chisq.test(unstack(tbl, form = freq.positive ~ group)[1:4,]) ### no diff in proportion of positive coeff
# 
# t.test(freq.impact/nb.sp ~ group, data=tbl, paired = T)  ### signif P= 0.01681, df=4, t=-3.9508
# t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T)   ### signif P = 0.03048, df=4, t=-32808


## other method : Independence tests : NOTHING SIGNIFICANT
tbl <- as.table(t(rbind(none = as.numeric(table(na=is.na(glmSRnat$boot.thresh$th), alien = rownames(glmSRnat$boot.thresh)%in%aliens)[2,]),
                        table(glmSRnat$boot.thresh$th, alien = rownames(glmSRnat$boot.thresh)%in%aliens))
))
library(coin)
cmh_test(tbl)
chisq_test(as.table( tbl))

