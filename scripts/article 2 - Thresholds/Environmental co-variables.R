### Testing effect of environmental factors 

focal.aliens <- rownames(glmSRnat.overall$glms)
focal.aliens <- focal.aliens[focal.aliens %in% aliens]

envir.covar <- data.frame(matrix(NA, 3 + length(focal.aliens), 14))
rownames(envir.covar) <- c("SR", "SRnat","SRali", focal.aliens)
colnames(envir.covar) <- c("signif","elev.R","elev.df","elev.P", 
                           "asp.R","asp.df","asp.P",
                           "elev.slope", "elev.P", 
                           "asp.slope", "asp.P",
                           "inter.slope", "inter.P",
                           "total.R2")

envir.covar$signif[rownames(envir.covar) %in% impsp] <- 1

## main drivers of Native Richness: ####
par(mfrow = c(4,5))
for (i in 2:20) {
  var <- names(envplot@data)[i]
  plot(envplot@data[,var],  envplot@data[,"SRnat"],  ann = F)
  fit <- cor.test(envplot@data[,var],  envplot@data[,"SRnat"],method = "pearson")
  envir.covar[i,1:3] <- c(fit$estimate, fit$parameter, fit$p.value)
  mtext(3, text = paste(round(fit$estimate, 2), p2star(fit$p.value)), adj = 1, cex = 0.7)
  mtext(3, text = var, adj = 0.5, cex = 0.7)
}


#### All species correlation with elevation ######

# elevation
par(mfrow = c(1,3))
for (i in 1:3) {
  var <- c("SR","SRnat", "SRali")[i] 
  plot(envplot$DEM_10, envplot@data[,var],  main = var, xlab = "elevation", ylab = "richness")
  fit <- cor.test(envplot$DEM_10, envplot@data[,var], method = "pearson")
  envir.covar[i,2:4] <- c(fit$estimate, fit$parameter, fit$p.value)
  mtext(3, text = paste(round(fit$estimate, 2), p2star(fit$p.value)), adj = 1)
}

for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- comm[, sp]
  abun[abun == 0] <- NA
  fit <- cor.test(envplot$DEM_10, abun, method = "pearson")
  envir.covar[i+3,2:4] <- c(fit$estimate, fit$parameter, fit$p.value)
}

# aspect
par(mfrow = c(1,3))
for (i in 1:3) {
  var <- c("SR","SRnat", "SRali")[i] 
  aspect <- abs(envplot$ASPECT - 180)
  plot(aspect, envplot@data[,var],  main = var, xlab = "Aspect", ylab = "species richness")
  fit <- cor.test(aspect, envplot@data[,var], method = "pearson")
  envir.covar[i,5:7] <- c(fit$estimate, fit$parameter, fit$p.value)
  mtext(3, text = paste(round(fit$estimate, 2), p2star(fit$p.value)), adj = 1)
}


par(mfrow = c(3,4))
for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- comm[, sp]
  abun[abun == 0] <- NA
  aspect <- abs(envplot$ASPECT - 180)
  fit <- cor.test(aspect, abun, method = "pearson")
  envir.covar[i +3,5:7] <- c(fit$estimate, fit$parameter, fit$p.value)
}



# aspect*elevation
par(mfcol = c(2,3))
for (i in 1:3) {
  var <- c("SR","SRnat", "SRali")[i] 
  aspect <- abs(envplot$ASPECT - 180) ## 0 = South, 180 = North
  altitude <- envplot$DEM_10
  plot(altitude, envplot@data[,var],  main = var, xlab = "Altitude", ylab = "species richness")
  plot(aspect, envplot@data[,var],  main = var, xlab = "Aspect", ylab = "species richness")
  fit1 <- lm(envplot@data[,var] ~ altitude)
  fit2 <- lm(envplot@data[,var] ~ altitude * aspect)
  if (AIC(fit2) < AIC(fit1) - 2 ) {
    fit <- fit2
    envir.covar[i,8:14] <- c(summary(fit)$coef[2,1], summary(fit)$coef[2,4], 
                           summary(fit)$coef[3,1], summary(fit)$coef[3,4], 
                           summary(fit)$coef[4,1], summary(fit)$coef[4,4], 
                           summary(fit)$adj.r.squared
                           )
  }
  else {
    fit <- fit1
    envir.covar[i,8:14] <- c(summary(fit)$coef[2,1], summary(fit)$coef[2,4], 
                             NA, NA,
                             NA, NA,
                             summary(fit)$adj.r.squared
                             )
  }
    
}


for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- comm[, sp]
  abun[abun == 0] <- NA
  aspect <- abs(envplot$ASPECT - 180)
  altitude <- envplot$DEM_10
plot(as.numeric(abun) ~ altitude)
  fit1 <- lm( abun ~ altitude)
  fit2 <- lm( abun ~ altitude * aspect)
  if (AIC(fit2) < AIC(fit1) - 2 ) {
    fit <- fit2
    envir.covar[i+3,8:14] <-c(summary(fit)$coef[2,1], summary(fit)$coef[2,4], 
                             summary(fit)$coef[3,1], summary(fit)$coef[3,4], 
                             summary(fit)$coef[4,1], summary(fit)$coef[4,4], 
                             summary(fit)$adj.r.squared
    )
  }
  else {
    fit <- fit1
    envir.covar[i+3,8:14] <- c(summary(fit)$coef[2,1], summary(fit)$coef[2,4], 
                             NA, NA,
                             NA, NA,
                             summary(fit)$adj.r.squared
    )
  }
  
}


envir.covar$elev.R2 <- envir.covar$elev.R^2
write.csv(envir.covar, file = "envir.covar.csv")

# number of species increasing with elevation: 6
sum(envir.covar[4:nrow(envir.covar),]$elev.R > 0.1 & envir.covar[4:nrow(envir.covar),]$elev.P<0.05, na.rm = T)
# number of species decreasing with elevation: 10
sum(envir.covar[4:nrow(envir.covar),]$elev.R < (-0.1) & envir.covar[4:nrow(envir.covar),]$elev.P<0.05, na.rm = T)

# number of species with significant interactions:
envir.covar [which(envir.covar$inter.P < 0.05 & envir.covar$total.R2 >0.01),]



### add spatial structure ####

# geographical distance matrix
library(raster)
spDst <-pointDistance(coordinates(envplot), lonlat = F, allpairs=FALSE)
rownames(spDst) <- colnames(spDst) <-rownames(coordinates(envplot))



###### Testing ABUN s a factor + elevation and aspect as a covariable for the glm #####
glms.covar <- data.frame(matrix(NA, length(focal.aliens), 12),row.names = focal.aliens)
colnames(glms.covar) <- c("df", "P.abun.alone","DeltaAIC", "P.elev", "P.asp", "P.abun", "est.abun", "P.c2", "P.c3", "P.c4", "P.c5", "P.c6")

# par(mfcol = c(1,2), mar = c(2,2,1,1), oma= c(1,1,3,1))
for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- comm[, sp]
  abun[abun == 0] <- NA
  elev <- envplot$DEM_10[!is.na(abun)]
  aspect <- abs(envplot$ASPECT - 180)[!is.na(abun)]
  Y <- envplot$SRnat[!is.na(abun)]
  abun <- na.omit(abun)
  # plot(elev, abun, ann = F)
  # mtext(3, text = sp)
  # plot(abun, Y, ann = F)
  fit0 <- glm( Y ~ elev + aspect, family = poisson)
  fit1 <- glm( Y ~ elev + aspect + as.factor(abun), family = poisson)
  fit2 <- glm( Y ~ elev + aspect + as.numeric(abun), family = poisson)
  fit3 <- glm( Y ~ as.numeric(abun), family = poisson)
  glms.covar[i,1:(7 + range(abun)[2]-1)] <- c(fit2$df.residual, 
                                              summary(fit3)$coef[2,4],
                                              AIC(fit0) - AIC(fit2), 
                                              summary(fit2)$coeff[2:4,4],
                                              summary(fit2)$coeff[4,1],
                                              summary(fit1)$coeff[4:length(coef(fit1)), 4])
  
  mtext(round(AIC(fit0) - AIC(fit1), 1), 3, adj = 0)
  mtext( paste("elev", p2star(summary(fit2)$coeff[2,4]),
               "asp", p2star(summary(fit2)$coeff[3,4]),
              "abun", p2star(summary(fit2)$coeff[4,4])))
  mtext(3, at = 1:5, line= -1, text = sapply(summary(fit1)$coeff[4:nrow(summary(fit1)$coeff), 4], FUN = p2star))
  lines(1:5, y = coef(fit2)[4] * (1:5) + coef(fit2)[1], col = "red")
} # RYTRAC and CRIMUR become non significant when adding altitude.



## Plotting the new 10 significant species : #####
quartz()
par(mfrow = c(5,4), mar = c(3,3,3,1), oma= c(1,1,1,1), las = 1, mgp = c(1.3,0.5,0), cex.axis = 0.8)
for (i in 1:length(impsp)) {
  sp <- impsp[i] 
  abun <- comm[, sp]
  abun[abun == 0] <- NA
  elev <- envplot$DEM_10[!is.na(abun)]
  aspect <- abs(envplot$ASPECT - 180)[!is.na(abun)]
  SRnat <- envplot$SRnat[!is.na(abun)]
  abun <- na.omit(abun)
  plot(elev, abun)
  mtext(3, text = sp, line= 1)
  plot(abun, SRnat)
  fit0 <- glm( SRnat ~ elev + aspect, family = poisson)
  fit1 <- glm( SRnat ~ elev + aspect + as.factor(abun), family = poisson)
  fit2 <- glm( SRnat ~ elev + aspect + as.numeric(abun), family = poisson)
  fit3 <- glm( SRnat ~ as.numeric(abun), family = poisson)
  mtext(paste("DeltaAIC:",round(AIC(fit0) - AIC(fit1), 1)), 3, adj = 0, cex = 0.7, line= 1.2)
  mtext( paste("elev", p2star(summary(fit2)$coeff[2,4]),
               "asp", p2star(summary(fit2)$coeff[3,4]),
               "abun", p2star(summary(fit2)$coeff[4,4])), cex = 0.7, adj = 0)
  p <- summary(fit1)$coeff[4:nrow(summary(fit1)$coeff), 4]
  p[summary(fit1)$coeff[4:nrow(summary(fit1)$coeff), 1] >0] <- NA
  
  mtext(3, at = 1:5, line= -2, text = sapply(p, FUN = p2star), col = "red")
  lines(1:5, y = coef(fit2)[4] * (1:5) + coef(fit2)[1], col = "red")
} 
 
 
 
# significance and with bootstrapping:
(sign(glmSRnat.overall$CIhi)*sign(glmSRnat.overall$CIlow) == 1)[focal.aliens,]


  
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

