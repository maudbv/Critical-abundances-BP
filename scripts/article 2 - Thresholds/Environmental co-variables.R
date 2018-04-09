### Testing effect of environmental factors 
focal.aliens <- rownames(glmSRnat.overall$glms)
focal.aliens <- focal.aliens[focal.aliens %in% aliens]

tmp.comm <- comm[unimprovedgrasslands,]
tmp.envplot <- envplot[unimprovedgrasslands,]


# map environmental variables?


#### CHoose a subset of environmental variables for analyses:
tmp=envplot@data[,c(  "PANN0080","PSUM0080", "PWIN0080",  "AVTEMP0080","MINTEM0080", "MAXTEM0080","GDD","SOLRADR","MOISTR","UR1991_DNS","BLDG_DIST","DIST_HWSLD","DIST_METAL","DIST_PRIV","DIST_SRIV", "DEM_10","ASPECT","SLOPE","PH_MID")]
tmp$dist_built = apply(tmp[, c("BLDG_DIST","DIST_HWSLD")], 1, min)
# tmp=tmp-colMeans(tmp)/ apply(tmp, 2, sd, na.rm=T)    
tmp=cbind(tmp,envplot@data[,c("Northern", "Eastern", "NWestern","SR","SRnat","SRali", "vegtype")])
tmp$logSRnat <- log(tmp$SRnat)
tmp$logSRnat[which(is.infinite(tmp$logSRnat))] <- NA


## PCA on environmental variables #####

library(FactoMineR)

                          
## PCA on envir variables
bb=PCA(tmp[,c( "vegtype","SRnat","SRali",
               "PANN0080", "AVTEMP0080","SOLRADR","MOISTR","DEM_10",
               "PH_MID","BLDG_DIST","DIST_HWSLD","Northern", "Eastern", "SLOPE")],
       quanti.sup=c(2,3), quali.sup=1)

b=PCA(tmp[,c( "SRnat" ,"SRali", "PANN0080", "AVTEMP0080","DEM_10",
              "Northern","BLDG_DIST","DIST_HWSLD", "SLOPE")],quanti.sup=c(1:4))
quartz()
par(mfrow=c(1,2))
plot(b, axes = c(1,2), choix= "var")
plot(b, axes = c(1,3), choix= "var")

plot( Northern ~ ASPECT, tmp)
plot( Eastern ~ ASPECT, tmp)
plot( Eastern ~ Northern, tmp)

plot( SRali  ~ DIST_HWSLD, tmp)
plot( SRnat  ~ DIST_HWSLD, tmp)
plot( SRnat  ~ BLDG_DIST, tmp)


plot( BLDG_DIST  ~ DIST_HWSLD, tmp)   ### Very correlated
cor.test(tmp$BLDG_DIST, tmp$DIST_HWSLD)

plot( SRnat  ~ DIST_HWSLD, tmp)   ### R2 = 12%
cor.test(tmp$SRnat, tmp$DIST_HWSLD)

plot( SRnat  ~ BLDG_DIST , tmp)   ### R2 = 19%
cor.test(tmp$SRnat, tmp$BLDG_DIST )

plot( logSRnat ~ dist_built, tmp)   ### Very correlated



## Variable selection 
library(leaps)

leaps<-regsubsets( logSRnat ~ DEM_10 + SLOPE + Northern + dist_built + SRali, data = tmp, nbest=10)
# view results
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
# plot statistic by subset size
library(car)
subsets(leaps, statistic="rsq") 




## Testing interaction of envir factors on SRnat =  GLM sith interactions ####
  tmp =  envplot@data[unimprovedgrasslands,]
 f0<- glm( SRnat ~ 1, data = tmp, family = poisson)
 f1<- glm( SRnat ~ DEM_10, data = tmp, family = poisson)
 f2<- glm( SRnat ~ DEM_10 * SLOPE, data = tmp, family = poisson)
 f3<- glm( SRnat ~ DEM_10 * SLOPE * Northern, data = tmp, family = poisson)
 f3ni<- glm( SRnat ~ DEM_10 + SLOPE + Northern, data = tmp, family = poisson)
 f3_ExN<- glm( SRnat ~ DEM_10 + SLOPE + DEM_10:Northern + DEM_10:SLOPE, data = tmp, family = poisson)
 f4<- glm( SRnat ~ DEM_10 * SLOPE * Northern * SRali, data = tmp, family = poisson)
 f4ni<- glm( SRnat ~ DEM_10 +  SLOPE   + Northern+ SRali, data = tmp, family = poisson)
 f4_simp<- glm( SRnat ~  DEM_10 + SLOPE + SRali + DEM_10:Northern + DEM_10:SLOPE + DEM_10:SRali, data = tmp, family = poisson)
 f4_supersimp<- glm( SRnat ~  DEM_10 + SLOPE + SRali + Northern + DEM_10:SRali, data = tmp, family = poisson)
  
 
AIC(f0,f1,f2,f3,f3ni,f3_ExN, f4,f4ni, f4_simp, f4_supersimp)
 summary(f3) 
 summary(f3ni) 
 summary(f3_ExN) # best model has interactions between elevation and the two others.
 anova(f3)
 summary(f4) #best model including SRali has all interactions
 summary(f4ni) 
 summary(f4_simp) 
 summary(f4_supersimp)
 anova(f4_supersimp)
 
 # no interactions
 f0<- glm( SRnat ~ 1, data = tmp, family = poisson)
 f1<- glm( SRnat ~ DEM_10, data = tmp, family = poisson)
 f2<- glm( SRnat ~ DEM_10 + SLOPE, data = tmp, family = poisson)
 f3<- glm( SRnat ~ DEM_10 + SLOPE + Northern, data = tmp, family = poisson)
f4<- glm( SRnat ~ DEM_10 + SLOPE + Northern + DIST_HWSLD , data = tmp, family = poisson)
f5<- glm( SRnat ~ DEM_10 + SLOPE + Northern + DIST_HWSLD +  BLDG_DIST, data = tmp, family = poisson)
f6<- glm( SRnat ~ DEM_10 + SLOPE + Northern + DIST_HWSLD +  BLDG_DIST + SRali, data = tmp, family = poisson)

f4simp<- glm( SRnat ~ DEM_10 + SLOPE + Northern + SRali , data = tmp, family = poisson)
anova(f4simp)


AIC(f0,f1,f2,f3,f4,f4simp, f5,f6)
anova(f3,f4,f5,f6)



 #For SRali
 f0<- glm( SRali ~ 1, data = tmp, family = poisson)
 f1<- glm( SRali ~ DEM_10, data = tmp, family = poisson)
 f2<- glm( SRali ~ DEM_10 * SLOPE, data = tmp, family = poisson)
 f3ni <-  glm( SRali ~ DEM_10 + SLOPE + Northern, data = tmp, family = poisson)
 f3<- glm( SRali ~ DEM_10 * SLOPE * Northern, data = tmp, family = poisson)
 
 AIC(f0,f1,f2,f3ni,f3) 
 anova(f3ni)
 summary(f3ni) # best model has no interactions
 
 
 # not interactions
 f0<- glm( SRali ~ 1, data = tmp, family = poisson)
 f1<- glm( SRali ~ DEM_10, data = tmp, family = poisson)
 f2<- glm( SRali ~ DEM_10 + SLOPE, data = tmp, family = poisson)
 f3<- glm( SRali ~ DEM_10 + SLOPE + Northern, data = tmp, family = poisson)
 f3<- glm( SRali ~ DEM_10 + SLOPE + Northern + SRnat, data = tmp, family = poisson)
 
AIC(f0,f1,f2,f3, f4) 
anova(f3)
summary(f3) # best model has no interactions




###  Creating a table ENVIR.COVAR summarizing all linear correlations between factors: ####

envir.covar <- data.frame(matrix(NA, 3 + length(focal.aliens), 20))
rownames(envir.covar) <- c("SR", "SRnat","SRali", focal.aliens)
colnames(envir.covar) <- c("signif","elev.R","elev.df","elev.P", 
                           "slope.R","slope.df","slope.P",
                           "north.R","north.df","north.P",
                           "SRali.R","SRali.df","SRali.P",
                           "elev.slope", "elev.slope.P", 
                           "north.slope", "north.slope.P",
                           "inter.slope", "inter.P",
                           "total.R2")

envir.covar$signif[rownames(envir.covar) %in% impsp] <- 1

# Elevation
par(mfrow = c(1,3))
for (i in 1:3) {
  var <- c("SR","SRnat", "SRali")[i] 
  plot(tmp.envplot$DEM_10, tmp.envplot@data[,var],  main = var, xlab = "elevation", ylab = "richness")
  fit <- cor.test(tmp.envplot$DEM_10, tmp.envplot@data[,var], method = "pearson")
  envir.covar[i,c("elev.R","elev.df","elev.P")] <- c(fit$estimate, fit$parameter, fit$p.value)
  mtext(3, text = paste(round(fit$estimate, 2), p2star(fit$p.value)), adj = 1)
}

for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- tmp.comm[, sp]
  abun[abun == 0] <- NA
  fit <- cor.test(tmp.envplot$DEM_10, abun, method = "pearson")
  envir.covar[i+3,c("elev.R","elev.df","elev.P")] <- c(fit$estimate, fit$parameter, fit$p.value)
}


# Slope
par(mfrow = c(1,3))
for (i in 1:3) {
  var <- c("SR","SRnat", "SRali")[i] 
  plot(tmp.envplot$SLOPE, tmp.envplot@data[,var],  main = var, xlab = "elevation", ylab = "richness")
  fit <- cor.test(tmp.envplot$SLOPE, tmp.envplot@data[,var], method = "pearson")
  envir.covar[i,c("slope.R","slope.df","slope.P")] <- c(fit$estimate, fit$parameter, fit$p.value)
  mtext(3, text = paste(round(fit$estimate, 2), p2star(fit$p.value)), adj = 1)
}

for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- tmp.comm[, sp]
  abun[abun == 0] <- NA
  fit <- cor.test(tmp.envplot$SLOPE, abun, method = "pearson")
  envir.covar[i+3,c("slope.R","slope.df","slope.P")] <- c(fit$estimate, fit$parameter, fit$p.value)
}

# Northern aspect
par(mfrow = c(1,3))
for (i in 1:3) {
  var <- c("SR","SRnat", "SRali")[i] 
  north <- tmp.envplot$Northern
  plot(north, tmp.envplot@data[,var],  main = var, xlab = "Northernness", ylab = "species richness")
  fit <- cor.test(north, tmp.envplot@data[,var], method = "pearson")
  envir.covar[i,c("north.R","north.df","north.P")] <- c(fit$estimate, fit$parameter, fit$p.value)
  mtext(3, text = paste(round(fit$estimate, 2), p2star(fit$p.value)), adj = 1)
}


par(mfrow = c(3,4))
for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- tmp.comm[, sp]
  abun[abun == 0] <- NA
  north <- tmp.envplot$Northern
  fit <- cor.test(north, abun, method = "pearson")
  envir.covar[i +3,c("north.R","north.df","north.P")] <- c(fit$estimate, fit$parameter, fit$p.value)
}



# SRali
par(mfrow = c(1,3))
for (i in 1:3) {
  var <- c("SR","SRnat", "SRali")[i] 
  plot(tmp.envplot@data[,"SRali"], tmp.envplot@data[,var],  main = var, xlab = "SRali", ylab = "species richness")
  fit <- cor.test(tmp.envplot@data[,"SRali"], tmp.envplot@data[,var], method = "pearson")
  envir.covar[i,c("SRali.R","SRali.df","SRali.P")] <- c(fit$estimate, fit$parameter, fit$p.value)
  mtext(3, text = paste(round(fit$estimate, 2), p2star(fit$p.value)), adj = 1)
}


par(mfrow = c(2,4))
for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- tmp.comm[, sp]
  abun[abun == 0] <- NA
  SRali <- tmp.envplot@data[,"SRali"]
  
  if (sp %in% impsp) {
    plot( abun,SRali,  main = sp, xlab = "SRali", ylab = "Abun")
  }
  fit <- cor.test(SRali, abun, method = "pearson")
  envir.covar[i +3,c("SRali.R","SRali.df","SRali.P")] <- c(fit$estimate, fit$parameter, fit$p.value)
}

# Northern*elevation
par(mfcol = c(2,3))
for (i in 1:3) {
  var <- c("SR","SRnat", "SRali")[i] 
  north <- tmp.envplot$Northern
  altitude <- tmp.envplot$DEM_10
  plot(altitude, tmp.envplot@data[,var],  main = var, xlab = "Altitude", ylab = "species richness")
  plot(north, tmp.envplot@data[,var],  main = var, xlab = "Northern", ylab = "species richness")
  fit1 <- lm(tmp.envplot@data[,var] ~ altitude)
  fit2 <- lm(tmp.envplot@data[,var] ~ altitude * north)
  if (AIC(fit2) < (AIC(fit1) - 2) ) {
    fit <- fit2
    envir.covar[i,c("elev.slope", "elev.slope.P","north.slope", "north.slope.P","inter.slope", "inter.P","total.R2")] <- c(summary(fit)$coef[2,1], summary(fit)$coef[2,4], 
                                                                                                                           summary(fit)$coef[3,1], summary(fit)$coef[3,4], 
                                                                                                                           summary(fit)$coef[4,1], summary(fit)$coef[4,4], 
                                                                                                                           summary(fit)$adj.r.squared
    )
  }
  else {
    fit <- fit1
    envir.covar[i,c("elev.slope", "elev.slope.P","north.slope", "north.slope.P","inter.slope", "inter.P","total.R2")] <- c(summary(fit)$coef[2,1], summary(fit)$coef[2,4], 
                                                                                                                           NA, NA,
                                                                                                                           NA, NA,
                                                                                                                           summary(fit)$adj.r.squared
    )
  }
  
}


for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- tmp.comm[, sp]
  abun[abun == 0] <- NA
  north <- tmp.envplot$Northern
  altitude <- tmp.envplot$DEM_10
  plot(as.numeric(abun) ~ altitude)
  fit1 <- lm( abun ~ altitude)
  fit2 <- lm( abun ~ altitude * north)
  if (AIC(fit2) < AIC(fit1) - 2 ) {
    fit <- fit2
    envir.covar[i+3,c("elev.slope", "elev.slope.P","north.slope", "north.slope.P","inter.slope", "inter.P","total.R2")] <-c(summary(fit)$coef[2,1], summary(fit)$coef[2,4], 
                                                                                                                            summary(fit)$coef[3,1], summary(fit)$coef[3,4], 
                                                                                                                            summary(fit)$coef[4,1], summary(fit)$coef[4,4], 
                                                                                                                            summary(fit)$adj.r.squared
    )
  }
  else {
    fit <- fit1
    envir.covar[i+3,c("elev.slope", "elev.slope.P","north.slope", "north.slope.P","inter.slope", "inter.P","total.R2")] <- c(summary(fit)$coef[2,1], summary(fit)$coef[2,4], 
                                                                                                                             NA, NA,
                                                                                                                             NA, NA,
                                                                                                                             summary(fit)$adj.r.squared
    )
  }
  
}


envir.covar$elev.R2 <- envir.covar$elev.R^2
write.csv(envir.covar, file = "envir.covar.csv")

# number of species increasing with elevation: 6
sum(envir.covar[4:nrow(envir.covar),]$elev.R > 0 & envir.covar[4:nrow(envir.covar),]$elev.P<=0.05, na.rm = T)
# number of species decreasing with elevation: 6
sum(envir.covar[4:nrow(envir.covar),]$elev.R < 0 & envir.covar[4:nrow(envir.covar),]$elev.P<=0.05, na.rm = T)


# number of species with some correlation to one of the gradients:
sum(rowSums(
  cbind(
    (envir.covar[4:nrow(envir.covar),]$elev.P<=0.05),
    (envir.covar[4:nrow(envir.covar),]$slope.P<=0.05),
    (envir.covar[4:nrow(envir.covar),]$north.P<=0.05)
  ),
  na.rm = T)>0)
23/47 # number of species with a trend along one of the environemntal gradients (not including SRali) = roughly half

# number of species with focal abundance opposite to native richenss patterns in relation to environmental gradients: 14
sum(rowSums(
  cbind(
    (envir.covar[4:nrow(envir.covar),]$elev.R <0 & envir.covar[4:nrow(envir.covar),]$elev.P<=0.05), # decrease with elevation
    (envir.covar[4:nrow(envir.covar),]$slope.R < 0 & envir.covar[4:nrow(envir.covar),]$slope.P<=0.05), # decrease with slope
    (envir.covar[4:nrow(envir.covar),]$north.R > 0 & envir.covar[4:nrow(envir.covar),]$north.P<=0.05) # increase with northness
  ),
    na.rm = T)>0)

# number of species with focal abundance similar to native richenss patterns in relation to environmental gradients:11
sum(rowSums(
  cbind(
    (envir.covar[4:nrow(envir.covar),]$elev.R >0 & envir.covar[4:nrow(envir.covar),]$elev.P<=0.05), # increases with elevation
    (envir.covar[4:nrow(envir.covar),]$slope.R > 0 & envir.covar[4:nrow(envir.covar),]$slope.P<=0.05),# increases with slope
    (envir.covar[4:nrow(envir.covar),]$north.R < 0 & envir.covar[4:nrow(envir.covar),]$north.P<=0.05) # decreases with northness
  ),
  na.rm = T)>0)

# number of species with significant interactions:
envir.covar [which(envir.covar$inter.P < 0.05 & envir.covar$total.R2 >0.01),]









### TABLES1 for focal species abundance correlation with covariables #####
tableS1.focalsp <- data.frame(matrix(NA, length(focal.aliens), 17),row.names = focal.aliens)
colnames(tableS1.focalsp) <- c("df","SS","RSS","P", "R2",
                               "elev.est","slope.est","north.est","SRali.est",
                               "elev.R2","slope.R2","north.R2","SRali.R2",
                               "elev.P","slope.P","north.P","SRali.P")


for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- tmp.comm[, sp]
  abun[abun == 0] <- NA
  elev <- tmp.envplot$DEM_10[!is.na(abun)]
  slope <- tmp.envplot$SLOPE[!is.na(abun)]
  SRali <- tmp.envplot$SRali[!is.na(abun)]
  SRnat <- tmp.envplot$SRnat[!is.na(abun)]
  north <- tmp.envplot$Northern[!is.na(abun)]
  abun <- na.omit(abun)

  fit0 <- lm( abun ~ 1) 
  fit4 <- lm( as.numeric(abun) ~   north  + slope + SRali + elev) 
  

  tableS1.focalsp[i, ] <- c(summary(fit4)$df[2],as.numeric(anova(fit0, fit4)[2,c(2,4,6)]),summary(fit4)$adj,
                            summary(fit4)$coef[2:5,1],
                            anova(fit4)[1:4,2]/sum(anova(fit4)[,2]), 
                            anova(fit4)[1:4,5])
}
 

write.csv(tableS1.focalsp, file= "table S1 for focal species.csv")

### TABLES1 WITH LOGIT Function for focal species abundance correlation with covariables #####

library(MASS)
tableS1.focalsp <- data.frame(matrix(NA, length(focal.aliens), 14),row.names = focal.aliens)
colnames(tableS1.focalsp) <- c("species", "resid.df","%dev.best","DeltaAIC.best", "P.best", "best",
                               "elev.OR","slope.OR","north.OR","SRali.OR",
                               "elev.P","slope.P","north.P","SRali.P")

options(contrasts = c("contr.treatment", "contr.poly"))

for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- tmp.comm[, sp]
  abun[abun == 0] <- NA
  elev <- tmp.envplot$DEM_10[!is.na(abun)]
  slope <- tmp.envplot$SLOPE[!is.na(abun)]
  SRali <- tmp.envplot$SRali[!is.na(abun)]
  SRnat <- tmp.envplot$SRnat[!is.na(abun)]
  north <- tmp.envplot$Northern[!is.na(abun)]
  abun <- na.omit(abun)

  fit0 <- polr(as.factor(abun) ~ 1,  Hess=TRUE)
  fit4 <- polr(as.factor(abun) ~elev  + slope +  north + SRali ,  Hess=TRUE)
  fit4.2 <- stepAIC(fit4, ~.^2, scope = list(upper = ~ elev  + slope +  north + SRali, lower = ~1))
  AOV <- anova(fit0,fit4.2, fit4)
 n.var = length(summary(fit4.2)$term)
 if (n.var = 4) AOV = anova(fit0, fit4)
 #coeffs
 (ctable <- coef(summary(fit4)))
 #Odds ratio = exp(coef) for a logit function
 OR <- exp(ctable)[1:4, 1]
 ## calculate and store p values
 p <- (pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2)[1:4]

 #output
 tableS1.focalsp[i, ] <- c(species[sp, "SpeciesName"],
                           fit4.2$df.residual,(AOV[1,3] - AOV[2,3])/AOV[1,3] ,AIC(fit0) - AIC(fit4.2),
                           AOV[2,7],paste(unlist(fit4.2$terms), collapse = " "),
                            OR,p)
}


write.csv(tableS1.focalsp, file= "table S1 for focal species.csv")


###### TABLES2 Testing ABUN s a factor + elevation+ SLOPE + aspect + SRALI  as a covariable for the glm #####
glms.covar <- data.frame(matrix(NA, length(focal.aliens), 18),row.names = focal.aliens)
colnames(glms.covar) <- c("df.alone", "percdev.abun.alone","deltaAIC.alone",
                          "df",
                          "est.elev","P.elev","percdev.elev",
                          "est.slope","P.slope","percdev.slope",
                          "est.north","P.north","percdev.north",
                          "est.SRali","P.SRali","percdev.SRali",
                          "percdev.abun",
                          "total.deltaAIC")

tmp.comm <- comm[unimprovedgrasslands,]
tmp.envplot <- envplot@data[unimprovedgrasslands,]
for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- tmp.comm[, sp]
  abun[abun == 0] <- NA
  elev <- tmp.envplot$DEM_10[!is.na(abun)]
  slope <- tmp.envplot$SLOPE[!is.na(abun)]
  SRali <- tmp.envplot$SRali[!is.na(abun)]
  SRnat <- tmp.envplot$SRnat[!is.na(abun)]
  north <- tmp.envplot$Northern[!is.na(abun)]
  abun <- na.omit(abun)
  # plot(elev, abun, ann = F)
  # mtext(3, text = sp)
  # plot(abun, Y, ann = F)
  fit0 <- glm( SRnat ~ 1 , family = poisson)
  # fit1 <- glm( SRnat ~ elev , family = poisson)
  # fit2 <- glm( SRnat ~ elev + slope, family = poisson) 
  # fit3 <- glm( SRnat ~ elev + slope + north, family = poisson)
   fit4 <- glm( SRnat ~ elev + slope + north + SRali, family = poisson) 
  fit5 <- glm( SRnat ~ elev + slope + north + SRali + as.factor(abun), family = poisson) 
  fit.alone <- glm( SRnat ~ as.factor(abun), family = poisson)
  
  aic.scores <- AIC(fit0,fit.alone,fit4,fit5)
  
  anova.fit5 <- anova(fit5)
  
  glms.covar[i,1:3] <- c(fit.alone$df.residual,
                          anova(fit.alone)[2,2]/anova(fit.alone)[1,4],
                         (aic.scores[1,2] - aic.scores[2,2]))
  
  glms.covar[i,4:18] <- c(anova.fit5[6,3],
                          summary(fit5)$coef[2,1],summary(fit5)$coef[2,4],
                          anova.fit5[2,2]/anova.fit5[1,4],
                          summary(fit5)$coef[3,1],summary(fit5)$coef[3,4],
                          anova.fit5[3,2]/anova.fit5[1,4],
                          summary(fit5)$coef[4,1],summary(fit5)$coef[4,4],
                          anova.fit5[4,2]/anova.fit5[1,4],
                          summary(fit5)$coef[5,1],summary(fit5)$coef[5,4],
                          anova.fit5[5,2]/anova.fit5[1,4],
                          anova.fit5[6,2]/anova.fit5[1,4],
                         (aic.scores[3,2] - aic.scores[4,2]))
  
} 

glms.covar$dev.ratio.actual <- tableS2[focal.aliens,]$dev.ratio
glms.covar$deltaAIC.actual <-tableS2[focal.aliens,]$aic.SRali- tableS2[focal.aliens,]$aic.abun

# add coefficient testing (negative and signif in bootstrap)
# significance and with bootstrapping:
total.tableS2 <- cbind(glms.covar, coef = glmSRnat.overall$est[focal.aliens,],
                  signif.boot = (sign(glmSRnat.overall$CIhi)*sign(glmSRnat.overall$CIlow) == 1)[focal.aliens,],
                  coef.nocovar = glmSRnat.overall.nocovar$est[focal.aliens,],
                  signif.boot.nocovar = (sign(glmSRnat.overall.nocovar$CIhi)*sign(glmSRnat.overall.nocovar$CIlow) == 1)[focal.aliens,]
)


write.csv(total.tableS2, "tableS2.withdeviance.csv")







## Explore main individual drivers of Native Richness: ####
par(mfrow = c(5,5), mar = c(2,2,2,0))
for (i in c(2:20, 49, 50)) {
  var <- names(envplot@data)[i]
  plot(envplot@data[,var],  envplot@data[,"SRnat"],  ann = F)
  fit <- cor.test(envplot@data[,var],  envplot@data[,"SRnat"],method = "pearson")
  envir.covar[i,1:3] <- c(fit$estimate, fit$parameter, fit$p.value)
  mtext(3, text = paste(round(fit$estimate, 2), p2star(fit$p.value)), adj = 1, cex = 0.7)
  mtext(3, text = var, adj = 0.5, cex = 0.7)
}



## Plotting the new 8 significant species : #####
quartz()

par(mfrow = c(2,4), mar = c(3,3,3,1), oma= c(1,1,1,1), las = 1, mgp = c(1.3,0.5,0), cex.axis = 0.8)
for (i in 1:length(impsp)) {
  sp <- impsp[i] 
  abun <- comm[, sp]
  abun[abun == 0] <- NA
  elev <- envplot$DEM_10[!is.na(abun)]
  aspect <- abs(envplot$ASPECT - 180)[!is.na(abun)]
  SRnat <- envplot$SRnat[!is.na(abun)]
  SRali <- envplot$SRali[!is.na(abun)]
  abun <- na.omit(abun)
  plot(elev, abun)
  mtext(3, text = sp, line= 1)
  plot(abun, SRnat)
  fit0 <- glm( SRnat ~ elev + aspect + SRali, family = poisson)
  fit1 <- glm( SRnat ~ elev + aspect  + SRali + as.factor(abun), family = poisson) ## abun as factor
  fit2 <- glm( SRnat ~ elev + aspect  + SRali + as.numeric(abun), family = poisson) ## abun as numeric
  
  mtext(paste("DeltaAIC:",round(AIC(fit0) - AIC(fit1), 1)), 3, adj = 0, cex = 0.7, line= 1.2)
  mtext( paste("elev", p2star(summary(fit2)$coeff[2,4]),
               "asp", p2star(summary(fit2)$coeff[3,4]),
               "SRali", p2star(summary(fit2)$coeff[4,4]),
               "abun", p2star(summary(fit2)$coeff[5,4])), cex = 0.7, adj = 0)
  p <- summary(fit1)$coeff[4:nrow(summary(fit1)$coeff), 4]
  p[summary(fit1)$coeff[4:nrow(summary(fit1)$coeff), 1] >0] <- NA
  
  mtext(3, at = 1:5, line= -2, text = sapply(p, FUN = p2star), col = "red")
  lines(1:5, y = coef(fit2)[4] * (1:5) + coef(fit2)[1], col = "red")
} 








### add spatial structure ####

# geographical distance matrix
library(raster)
spDst <-pointDistance(coordinates(envplot), lonlat = F, allpairs=FALSE)
rownames(spDst) <- colnames(spDst) <-rownames(coordinates(envplot))






## GLM by layers of altitude ######

###### Testing elevation as a covariable for the glm

# Elevation bands
hist(envplot$DEM_10)
library(Hmisc)
elev.bands <- cut2(envplot$DEM_10, g = 6, m = 50)
elev.bands <- cut2(envplot$DEM_10, cuts = c(100,200,300, 400)) ## cf Tomasetto et al. 2013

for (i in 1:10) {
  par(mfrow = c(5,2), mar = c(2,2,1,1), oma= c(1,1,3,1))
  sp <- impsp[i] 
  abun <- comm[, sp]
  abun[abun == 0] <- NA
  
  for (j in levels(elev.bands)) {
  select <- which(elev.bands== j)
  plot(envplot$DEM_10[select], abun[select], ann = F)
  plot(abun[select], envplot$SRnat[select], ann = F)
  # fit <- glm(envplot$SRnat[select] ~ abun[select], family = poisson)
  # fit1 <- glm(envplot$SRnat[select] ~ abun[select] * envplot$ASPECT[select] , family = poisson)
  # print(sp)
  # if (AIC(fit) - AIC(fit1) > 2) fit <- fit1  
  # pred <- predict.glm(fit, type = "response") # ????
  # y <- fit$coefficients[1] + (1:5) * fit$coefficients[2]
  # lines(1:5, y)
  
  x = as.numeric(abun[select])
  y = envplot$SRnat[select]
  fit <- lm(y~x)
  pred <- as.data.frame(predict(fit, newdata = data.frame(x = c(1:6)), interval = "confidence"))
  lines(pred$fit)
  lines(pred[,2], lty = "dotted")
  lines(pred[,3], lty = "dotted")
  
  mtext(3, text = paste("s =",round(summary(fit)$coefficients[2,1], 2),
                        "%D =", round(1 - fit$deviance/fit$null.deviance, 2),
                        p2star(summary(fit)$coefficients[2,4])),
        cex = 0.7, adj = 1)

  # print(summary(fit0))
  # print(summary(fit1))
} 
  mtext(3, text = sp, outer= TRUE)
}


#### Test gamma with elevation


# Elevation bands

library(Hmisc)
elev.bands <- cut2(envplot$DEM_10, g = 6, m = 50)
elev.bands.lo <- as.numeric(sapply(levels(elev.bands), function(x) substr(x, start = 2, stop = 4)))

# elev.bands <- cut2(envplot$DEM_10, cuts = c(100,200,300,400)) ## cf Tomasetto et al. 2013

gammas.nat <- sapply(levels(elev.bands), function(j) {
    select <- which(elev.bands== j)
   # gama richness in the band
      g <- sum(colSums(comm[select,colnames(comm) %in% natives]) >0, na.rm = T)
      return(g)
    })
gammas.ali <- sapply(levels(elev.bands), function(j) {
  select <- which(elev.bands== j)
  # gama richness in the band
  g <- sum(colSums(comm[select,colnames(comm) %in% aliens]) >0, na.rm = T)
  return(g)
})
gammas.tot <- sapply(levels(elev.bands), function(j) {
  select <- which(elev.bands== j)
  # gama richness in the band
  g <- sum(colSums(comm[select,]) >0, na.rm = T)
  return(g)
})

plot(elev.bands.lo, gammas.tot,type = "b", xlab = "Elevational bands (m)",ylab = "Gamma richness", 
     ylim=c(0, 400), las = 1,  xaxt = "n")
axis(1, at = elev.bands.lo, labels = levels(elev.bands), cex.axis = 0.7)
points(elev.bands.lo, gammas.nat, type = "b",col = "forestgreen")
points(elev.bands.lo, gammas.ali, type = "b",col = "firebrick")
mtext(3, at = elev.bands.lo, text = summary(elev.bands), line=-0.7, cex = 0.7, col = "darkgrey")
legend("bottomleft", legend = c("All species", "Natives", "Aliens"), lty = c("solid","solid","solid"), col = c("black","forestgreen","firebrick"), bty ="n", cex = 0.8)


# Elevation bands
hist(envplot$DEM_10)
library(Hmisc)
elev.bands <- cut2(envplot$DEM_10, g = 6, m = 50)
elev.bands <- cut2(envplot$DEM_10, cuts = c(100,200,300, 400)) ## cf Tomasetto et al. 2013

par(mfrow = c(5,2), mar = c(2,2,1,1), oma= c(1,1,3,1))
cols = colorRampPalette(c("lightblue","firebrick")) (6)
for (i in 1:10) {
  sp <- impsp[i] 
plot(c(1,6), c(0,300), type = "n") 
  for (j in 1:length(levels(elev.bands))) {
    select <- which(elev.bands== levels(elev.bands)[j])
    mat <- comm[select,]
    abun <- mat[, sp]
    abun[abun == 0] <- NA
    
    # gama richness for each abun level:
    gammas <- sapply(sort(unique(abun)), FUN = function(x) {
     g <- sum(colSums(mat[ which(abun == x),]) >0)
     names(g) <- x
     return(g)
    })
    points(gammas, type = "b", col = cols[j])
  } 
  mtext(3, text = sp, outer=FALSE)
}

