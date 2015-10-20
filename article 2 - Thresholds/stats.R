## Statistical analyses


# summarize thresholds

# create overall contingency tables
threshold = "th.CI"
glm.table <- lapply(list(Total.richness = glmSR.overall, Native.richness = glmSRnat.overall, Alien.richness =glmSRali.overall), FUN= function (M) {
  tbl <- as.table(t(rbind(none = as.numeric(table(na=is.na(M$impact.spread[, threshold]), alien = rownames(M$impact.spread)%in%aliens)[2,]),
                          table(M$impact.spread[, threshold], alien = rownames(M$impact.spread)%in%aliens))))
  return(tbl)
})

threshold = "pth.CI"
glm.table.positives <- lapply(list(Total.richness = glmSR.overall, Native.richness = glmSRnat.overall, Alien.richness =glmSRali.overall), FUN= function (M) {
  tbl <- as.table(t(rbind(none = as.numeric(table(na=is.na(M$impact.spread[, threshold]), alien = rownames(M$impact.spread)%in%aliens)[2,]),
                          table(M$impact.spread[, threshold], alien = rownames(M$impact.spread)%in%aliens))))
  return(tbl)
})

## compare natives and aliens total threshold/impact occurrence :

tbl<- glmSR.sum$overall.summary
tbl$nothr <- tbl$nb.sp - tbl$freq.thr
tbl= t(tbl[,3:4])
fisher.test(tbl)
chisq.test(tbl)

tbl <- glmSRnat.sum$overall.summary
tbl$nothr <- tbl$nb.sp - tbl$freq.thr
tbl= t(tbl[,3:4])
fisher.test(tbl)  #### the only interesting one -> not signif if more stringent on robustness...
chisq.test(tbl)

tbl <- glmSRali.sum$overall.summary
tbl$nothr <- tbl$nb.sp - tbl$freq.thr
tbl <- t(tbl[,3:4])
fisher.test(tbl)
chisq.test(tbl)


### MEdians
alien.ind <- which(rownames(glmSRnat.overall$impact.spread) %in% aliens)
median( glmSRnat.overall$impact.spread$th.CI[alien.ind], na.rm=T)

median( glmSRali.overall$impact.spread$th.CI[alien.ind], na.rm=T)



### better table for AR and NR
ali.tab <- data.frame(
  type = c(rep("NR",5), rep("AR",5)),
  class = c(2:6, 2:6),
  total.obs = c(glmSRnat.sum$class.summary$nb.sp[6:10],glmSRnat.sum$class.summary$nb.sp[6:10]),
  th = c(glmSRnat.sum$class.summary$freq.th[6:10],glmSRali.sum$class.summary$freq.th[6:10])
  )


# VGLM using VGAM package
library(VGAM)
# (thresh, no thresh) ~ class.abun + alien.status
#frequency of negative ES :
ali.tab <- data.frame(
  type = c(rep("NR",5), rep("AR",5)),
  class = c(2:6, 2:6),
  total.obs = c(glmSRnat.sum$class.summary$nb.sp[6:10],glmSRnat.sum$class.summary$nb.sp[6:10]),
  th = c(glmSRnat.sum$class.summary$freq.negative[6:10],glmSRali.sum$class.summary$freq.negative[6:10])
)

ali.tab$noth <- ali.tab$total.obs - ali.tab$th

vglm.SRnat<- vglm(cbind(th, noth) ~ class + type, data = ali.tab, family = binomialff,trace = TRUE)
summary(vglm.SR)


# (thresh, no thresh) ~ class.abun + type of richness
#frequency of negative ES :
ali.tab <- data.frame(
  type = c(rep("NR",5), rep("AR",5)),
  class = c(2:6, 2:6),
  total.obs = c(glmSRnat.sum$class.summary$nb.sp[6:10],glmSRnat.sum$class.summary$nb.sp[6:10]),
  th = c(glmSRnat.sum$class.summary$freq.negative[6:10],glmSRali.sum$class.summary$freq.negative[6:10])
)
vglm.SR<- vglm(cbind(th, noth) ~ class + type, data = ali.tab, family = binomialff,trace = TRUE)

# different proportion of negatives between alien and natives for SRnat
vglm.SRnat<- vglm(cbind(freq.positive, freq.negative) ~ class + group, data = glmSRnat.sum$class.summary, family = binomialff,trace = TRUE)
summary(vglm.SRnat)
  # => not significant for thereshold frequencies, probs because of zeros

# NO proportion of negatives between alien and natives for SRali
vglm.SRali<- vglm(cbind(freq.positive, freq.negative) ~ class + group, data = glmSRali.sum$class.summary, family = binomialff,trace = TRUE)
summary(vglm.SRali)


### counting number of species in common for three components of richness :
rownames(glmSRnat.overall$impact.spread)[which(!is.na(glmSRnat.overall$impact.spread$th) & !is.na(glmSR.overall$impact.spread$th))]
rownames(glmSRnat.overall$impact.spread)[which(is.na(glmSRnat.overall$impact.spread$th) & is.na(glmSR.overall$impact.spread$th))]
rownames(glmSRnat.overall$impact.spread)[which(!is.na(glmSRnat.overall$impact.spread$th) & !is.na(glmSRali.overall$impact.spread$th))]


## compare natives and aliens threshold occurrence per class:


#SRnat
tbl <- glmSRnat.sum $class.summary
t.test(freq.negative/nb.sp ~ group, data = tbl, paired = T) ### there is a difference in the proportion of negative coef per class
t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T)   ### signif P = 0.03048, df=4, t=-32808
t.test(freq.thr/freq.negative ~ group, data = tbl, paired =T)
round(tbl$freq.thr/tbl$nb.sp , 4)

tbl <- glmSRali.sum $class.summary
t.test(freq.negative/nb.sp ~ group, data = tbl, paired = T)
t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T)
t.test(freq.thr/freq.negative ~ group, data = tbl, paired =T)
round(tbl$freq.thr/tbl$nb.sp , 4)

## compare natives and aliens threshold occurrence per class: Mann whitney U two-sample test

tbl <- glmSR.sum $class.summary
wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$freq.th[tbl$group =="ALIEN:1"],
            alternative ="two.sided", paired=T)

tbl <- glmSRnat.sum $class.summary
print(wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$freq.th[tbl$group =="ALIEN:1"],
                  alternative ="two.sided", paired=T))


tbl <- glmSRali.sum $class.summary
wilcox.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$freq.th[tbl$group =="ALIEN:1"],
            alternative ="two.sided", paired=T)


# proportion test
tbl <- glmSRnat.sum$class.summary
prop.test(tbl$freq.negative[tbl$group =="ALIEN:1"], tbl$nb.sp[tbl$group =="ALIEN:1"])  # non random proportion of negative coeffs
prop.test(tbl$freq.negative[tbl$group =="ALIEN:0"], tbl$nb.sp[tbl$group =="ALIEN:0"]) # random for natives


prop.test(tbl$freq.th[tbl$group =="ALIEN:1"], tbl$nb.sp[tbl$group =="ALIEN:1"])  # random proportion of negative coeffs
prop.test(tbl$freq.negative.sig[tbl$group =="ALIEN:1"], tbl$nb.sp[tbl$group =="ALIEN:1"])  # non random proportion of negative coeffs
prop.test(tbl$freq.negative.above[tbl$group =="ALIEN:1"], tbl$nb.sp[tbl$group =="ALIEN:1"])  # non random proportion of negative coeffs
prop.test(tbl$freq.th[tbl$group =="ALIEN:1"], tbl$freq.negative[tbl$group =="ALIEN:1"])  # random proportion of negative coeffs


tbl <- glmSRali.sum$class.summary
prop.test(tbl$freq.negative[tbl$group =="ALIEN:1"], tbl$nb.sp[tbl$group =="ALIEN:1"])  # non random proportion of negative coeffs
prop.test(tbl$freq.negative[tbl$group =="ALIEN:0"], tbl$nb.sp[tbl$group =="ALIEN:0"]) # random for natives

prop.test(tbl$freq.th[tbl$group =="ALIEN:1"], tbl$nb.sp[tbl$group =="ALIEN:1"])  # non random proportion of negative coeffs
prop.test(tbl$freq.th[tbl$group =="ALIEN:0"], tbl$nb.sp[tbl$group =="ALIEN:0"]) # random for natives


## chisquare
chisq.test(unstack(tbl, form = freq.negative ~ group)) ### no difference in distrib of occurrence per class
t.test(n.target/nb.sp ~ group, data = tbl, paired = T) ### there is a difference in the proportion of target sp per class
t.test(freq.negative/nb.sp ~ group, data = tbl, paired = T) ### there is a difference in the proportion of target sp per class

chisq.test(unstack(tbl, form = freq.thr ~ group)[c(1,3,4),]) ### ns
chisq.test(unstack(tbl, form = freq.negative ~ group)) ### diff in proportion of negative coeff
chisq.test(unstack(tbl, form = freq.positive ~ group)[1:4,]) ### no diff in proportion of positive coeff

t.test(freq.impact/nb.sp ~ group, data=tbl, paired = T)  ### signif P= 0.01681, df=4, t=-3.9508
t.test(freq.thr/nb.sp ~ group, data = tbl, paired =T)   ### signif P = 0.03048, df=4, t=-32808

