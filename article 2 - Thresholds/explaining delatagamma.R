### Model explaining gamma richness with threshold + clustering

thresh <-  glmSRnat.overall$impact.spread[impsp,  "th.CI"]
cluster<- - (mnnd.th[impsp, "mnnd.obs_above"] - mnnd.th[impsp, "null.mean_above"])/ mnnd.th[impsp, "null.sd_above"] 
dg<- -(table.div.part$GRc - table.div.part$GRnull )

summary(f0 <- lm( dg ~ 1, constrast= "contr.helmert"))
summary(f1 <- lm( dg ~ cluster))
summary(f2 <- lm( dg ~ as.numeric(thresh)))
summary(f3 <- lm( dg ~ as.numeric(thresh) + cluster))

anova(f0,f1,  f3)
anova(f0,f2,  f3)
AIC(f0,f1,  f3)
AIC(f0,f2,  f3)
AIC(f1, f2)

