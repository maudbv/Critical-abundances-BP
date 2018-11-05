## TEst for bias in botanical surveys


library(data.table)
surveydata <- fread(file="data/survey metadata.csv")


surveydata$date <- as.Date(paste(surveydata$Day, surveydata$Month, surveydata$Year, sep = "-"), format = "%d-%m-%Y")

databp$date <- surveydata$date[match(databp$PlotFK, surveydata$PlotPK)]
databp$year <- as.numeric(format(databp$date, format = "%Y"))
databp$month <- as.numeric(format(databp$date, format = "%m"))
databp$month.t <- databp$month - 7 
databp$month.t[databp$month.t < (1)] <- databp$month [databp$month.t < (1)] + 5
plot(month.t ~month, databp)

databp$DOY <- as.numeric(format(databp$date, format = "%j"))
databp$DOY.t <- cos((databp$DOY * pi)/180)

envplot$date <- databp$date[match(as.character(envplot$PLOTID), as.character(databp$PlotName))]
envplot$year <- as.numeric(format(envplot$date, format = "%Y"))
envplot$month <- as.numeric(format(envplot$date, format = "%m"))
envplot$month.t <- envplot$month - 7 
envplot$month.t[envplot$month.t < (1)] <- envplot$month [envplot$month.t < (1)] + 5
plot(month.t ~month, envplot)

envplot$DOY <- as.numeric(format(envplot$date, format = "%j"))
envplot$DOY.t <- cos((envplot$DOY * pi)/180)

plot(envplot$DOY, envplot$DOY.t)


# transform for southern hemisphere
envplot$month.t <- abs(envplot$month -6)  # June = 0 ( dip of winter), December = 6 (peak of summer)


## Test correlation with SRnat

envplot$SR <- as.numeric(envplot$SR)

par(mfrow = c(1,3))
plot(SR ~ date, envplot, col = "grey")
summary(f <- glm(SR ~ as.numeric(year), envplot, family = poisson))
abline(lm(SR ~ date, envplot))
OR <- round(exp(coef(f)[2])*100 - 100, 2) 
mtext(3, text = paste( OR ,"% per year"))

plot(SRnat ~ date, envplot, col = "grey")
summary(f <- glm(SRnat ~ as.numeric(year), envplot, family = poisson))
abline(lm(SRnat ~ date, envplot))
OR <- round(exp(coef(f)[2])*100 - 100, 2) 
mtext(3, text = paste( OR ,"% per year"))

plot(SRali ~ date, envplot, col = "grey")
summary(f <- glm(SRali ~ as.numeric(year), envplot, family = poisson))
abline(lm(SRali ~ date, envplot))
OR <- round(exp(coef(f)[2])*100 - 100, 2) 
mtext(3, text = paste( OR ,"% per year"))


par(mfrow = c(1,3))
plot(SR ~ DOY.t, envplot, col = "grey", xaxt = "n", xlab = "cos(DOY)")
axis(1, at = c(-1,1), label = c( "Winter solstice", "Summer solstice"))
summary(f <- glm(SR ~ as.numeric(DOY.t), envplot, family = poisson))
f0 <-glm(SR ~ 1, envplot, family = poisson)
AIC(f0,f)
abline(lm(SR ~ DOY.t, envplot))
OR <- round(exp(coef(f)[2])*100 - 100, 2) 
mtext(3, text = paste( OR ,"% per DOY.t"))

plot(SRnat ~ DOY.t, envplot, col = "grey")
summary(f <- glm(SRnat ~ as.numeric(DOY.t), envplot, family = poisson))
OR <- round(exp(coef(f)[2])*100 - 100, 2) 
f0 <-glm(SRnat ~ 1, envplot, family = poisson)
AIC(f0,f)
abline(lm(SRnat ~ DOY.t, envplot))
mtext(3, text = paste( OR ,"% per DOY.t"))

plot(SRali ~ DOY.t, envplot, col = "grey")
summary(f <- glm(SRali ~ as.numeric(DOY.t), envplot, family = poisson))
f0 <-glm(SRali ~ 1, envplot, family = poisson)
AIC(f0,f)
abline(lm(SRali ~ DOY.t, envplot))
OR <- round(exp(coef(f)[2])*100 - 100, 2) 
mtext(3, text = paste( OR ,"% per DOY.t"))


##n SR nat per month
par(mfrow = c(1,3))
plot(SR ~ month.t, envplot, col = "grey", xaxt = "n", xlab = "cos(DOY)")
axis(1, at = c(-1,1), label = c( "Winter solstice", "Summer solstice"))
summary(f <- glm(SR ~ as.numeric(month.t), envplot, family = poisson))
f0 <-glm(SR ~ 1, envplot, family = poisson)
AIC(f0,f)
abline(lm(SR ~ month.t, envplot))
OR <- round(exp(coef(f)[2])*100 - 100, 2) 
mtext(3, text = paste( OR ,"% per month.t"))

plot(SRnat ~ month.t, envplot, col = "grey")
summary(f <- glm(SRnat ~ as.numeric(month.t), envplot, family = poisson))
OR <- round(exp(coef(f)[2])*100 - 100, 2) 
f0 <-glm(SRnat ~ 1, envplot, family = poisson)
AIC(f0,f)
abline(lm(SRnat ~ month.t, envplot))
mtext(3, text = paste( OR ,"% per month.t"))

plot(SRali ~ month.t, envplot, col = "grey")
summary(f <- glm(SRali ~ as.numeric(month.t), envplot, family = poisson))
f0 <-glm(SRali ~ 1, envplot, family = poisson)
AIC(f0,f)
abline(lm(SRali ~ month.t, envplot))
OR <- round(exp(coef(f)[2])*100 - 100, 2) 
mtext(3, text = paste( OR ,"% per month.t"))


## Test correlation with environment : no significant temporal pattern except SLOPE

envplot$SR <- as.numeric(envplot$SR)

par(mfrow = c(1,2))
plot(DEM_10 ~ date, envplot, col = "grey")
summary(f <- lm(DEM_10 ~ as.numeric(year), envplot))
plot(DEM_10 ~ DOY.t, envplot, col = "grey")
summary(f <- lm(DEM_10 ~ as.numeric(DOY.t), envplot))


par(mfrow = c(1,2))
plot(Northern ~ date, envplot, col = "grey")
summary(f <- lm(Northern~ as.numeric(year), envplot))
plot(Northern ~ DOY.t, envplot, col = "grey")
summary(f <- lm(Northern ~ as.numeric(DOY.t), envplot))


par(mfrow = c(1,2))
plot(SLOPE ~ date, envplot, col = "grey")  ### SIGNIF DECREASE IN SLOPE
summary(f <- lm(SLOPE~ as.numeric(year), envplot))
abline(lm(SLOPE ~ date, envplot))
plot(SLOPE ~ DOY.t, envplot, col = "grey")
summary(f <- lm(SLOPE ~ as.numeric(DOY.t), envplot))


## Test correlation with abundance of target species
source('~/Dropbox/Work/doc boulot/postdoc Berlin/R projects/utility functions/FUNCTION add stats.R', echo=TRUE)

par(mfrow = c(2,4), xpd = T, mar = c(3,1,1,1))
for (sp in impsp){
  dat <- databp[ databp$SpeciesCode == sp,]
  # barplot(table(dat$year,dat$abun) , beside = T )
  
  boxplot( abun ~ as.numeric(year), dat)
  add.stats(f <- lm(abun ~ as.numeric(year), dat))
  
  msr <- tapply(dat$SRnat,as.numeric(dat$year), mean, na.rm = T)
  ssr <-  tapply(dat$SRnat,as.numeric(dat$year), sd, na.rm = T)
  
  # plot(as.numeric(names(msr)), msr, ylim = c(-5,20), xlab = "months", ylab = "SRnat")
  # segments(as.numeric(names(msr)), msr+ssr,as.numeric(names(msr)), msr -ssr)
  # add.stats(f <- lm(SRnat ~ as.numeric(year), dat))
  # 
  # m <- tapply(dat$abun,as.numeric(dat$year), mean, na.rm = T)
  # s <-  tapply(dat$abun,as.numeric(dat$year), sd, na.rm = T)
  # n <- table(as.numeric(dat$year))
  
# 
#   plot(msr ~m, type = "n")
#   text(m, msr, labels = names(m))
#   add.stats(f <- lm(msr ~ m))
#   
}



# OVERALL Yearly means in richness and environmental variables explored
library(doBy)
y.envplot <- summaryBy(SR + SRnat + SRali + DEM_10 + SLOPE + Northern + Eastern ~ as.numeric(year),
                            data = as.data.frame(envplot)[envplot$PLOTID %in% unimprovedgrasslands,], FUN=c(mean,sd))
y.envplot$year<- as.numeric(y.envplot$year)     

par(mfrow = c(2,2), xpd = F, mar = c(4,4,2,2))

 plot( SR.mean ~ as.numeric(year), y.envplot, ylim = c(0,60), xlab = "year", ylab = "Species richness") 
 segments(as.numeric(y.envplot$year), y.envplot$SR.mean +  y.envplot$SR.sd, y.envplot$year, y.envplot$SR.mean -  y.envplot$SR.sd)
 points( SRnat.mean ~ as.numeric(year), y.envplot, col = "forestgreen") 
 segments(as.numeric(y.envplot$year), y.envplot$SRnat.mean +  y.envplot$SRnat.sd, y.envplot$year, y.envplot$SRnat.mean -  y.envplot$SRnat.sd,col = "forestgreen")
 points( SRali.mean ~ as.numeric(year), y.envplot,  col ="firebrick") 
 segments(as.numeric(y.envplot$year), y.envplot$SRali.mean +  y.envplot$SRali.sd, y.envplot$year, y.envplot$SRali.mean -  y.envplot$SRali.sd,col ="firebrick")
 add.stats(lm(SRnat.mean ~ as.numeric(year), y.envplot))
 legend("topright", legend = c("All species", "Alien", "Native"), col = c("black", "firebrick", "forestgreen"), pch = 1 )

 plot( SLOPE.mean ~ as.numeric(year), y.envplot,  ylim = c(0,40), xlab = "year", ylab = "Slope (degrees)") 
 segments(as.numeric(y.envplot$year), y.envplot$SLOPE.mean +  y.envplot$SLOPE.sd, y.envplot$year, y.envplot$SLOPE.mean -  y.envplot$SLOPE.sd)
 add.stats(lm(SLOPE.mean ~ as.numeric(year), y.envplot))
 
 plot( DEM_10.mean ~ as.numeric(year), y.envplot,  ylim = c(0,700), xlab = "year", ylab = "Elevation (m)") 
 segments(as.numeric(y.envplot$year), y.envplot$DEM_10.mean +  y.envplot$DEM_10.sd, y.envplot$year, y.envplot$DEM_10.mean -  y.envplot$DEM_10.sd)
 add.stats(lm(DEM_10.mean ~ as.numeric(year), y.envplot))
 
 plot( Northern.mean ~ as.numeric(year), y.envplot,  ylim = c(-1,1), xlab = "year", ylab = "Northness (cos(rad))") 
 segments(as.numeric(y.envplot$year), y.envplot$Northern.mean +  y.envplot$Northern.sd, y.envplot$year, y.envplot$Northern.mean -  y.envplot$Northern.sd)
 add.stats(lm(Northern.mean ~ as.numeric(year), y.envplot))

 ## Correlation between envir and SRnat over years
 par(mfrow = c(1,3), xpd = F, mar = c(4,4,2,2))
 
 plot( SRnat.mean ~ SLOPE.mean, y.envplot, ylim = c(0,20), xlab = "Mean Slope per year", ylab = "Mean Native Richness per year") 
 segments(as.numeric(y.envplot$SLOPE.mean), y.envplot$SRnat.mean +  y.envplot$SRnat.sd, y.envplot$SLOPE.mean, y.envplot$SRnat.mean -  y.envplot$SRnat.sd)
 add.stats(lm(SRnat.mean ~ SLOPE.mean, y.envplot))
 
 plot( SRnat.mean ~ DEM_10.mean, y.envplot, ylim = c(0,20), xlab = "Mean Elevation per year", ylab = "Mean Native Richness per year") 
 segments(as.numeric(y.envplot$DEM_10.mean), y.envplot$SRnat.mean +  y.envplot$SRnat.sd, y.envplot$DEM_10.mean, y.envplot$SRnat.mean -  y.envplot$SRnat.sd)
 add.stats(lm(SRnat.mean ~ DEM_10.mean, y.envplot))
 
 plot( SRnat.mean ~ Northern.mean, y.envplot, ylim = c(0,20), xlab = "Mean Northness per year", ylab = "Mean Native Richness per year") 
 segments(as.numeric(y.envplot$Northern.mean), y.envplot$SRnat.mean +  y.envplot$SRnat.sd, y.envplot$Northern.mean, y.envplot$SRnat.mean -  y.envplot$SRnat.sd)
 add.stats(lm(SRnat.mean ~ Northern.mean, y.envplot))
 
 
# Find best model : year is signifciant, but not very

library(bestglm)

f.bic <- bestglm( as.data.frame(envplot)[ unimprovedgrasslands, c("DOY.t","SLOPE","DEM_10", "Northern", "SRali","year", "SRnat")],family = poisson, IC = "BIC",RequireFullEnumerationQ = TRUE)
f.bic $BestModels
summary(f.bic$BestModel)
anova(f.bic$BestModel)

summary(f <- glm(SRnat ~ DEM_10 +SLOPE + Northern +year + SRali  , data = as.data.frame(envplot)[ unimprovedgrasslands,], family = poisson))  # Significant interaction of year and SRali: both richness increase together 
anova(f)

## INTERACTIONS??
summary(f <- glm(SRnat ~ DEM_10+SLOPE + year +  Northern +year:DEM_10  + DEM_10:SLOPE:year +year*SRali, data = as.data.frame(envplot)[ unimprovedgrasslands,], family = poisson))  # Significant interaction of year and SRali: both richness increase together 


# f.cv <- bestglm( as.data.frame(envplot)[ unimprovedgrasslands, c("year","month.t", "SLOPE","DEM_10", "Northern", "SRali","SRnat")],family = poisson, IC = "CV")

  
  
#investigate role of slope effect
par(mfrow = c(8,3), xpd = T, mar = c(3,3,2,1))
for (sp in impsp){
  dat <- databp[ databp$SpeciesCode == sp,]
  # barplot(table(dat$year,dat$abun) , beside = T )
  
  boxplot( SLOPE ~ as.numeric(year), dat)
  add.stats(f <- lm(abun ~ as.numeric(year), dat))
  
  msr <- tapply(dat$SRnat,as.numeric(dat$year), mean, na.rm = T)
  ssr <-  tapply(dat$SRnat,as.numeric(dat$year), sd, na.rm = T)
  
  plot(as.numeric(names(msr)), msr, ylim = c(-5,20), xlab = "months", ylab = "SRnat")
  segments(as.numeric(names(msr)), msr+ssr,as.numeric(names(msr)), msr -ssr)
  add.stats(f <- lm(SRnat ~ as.numeric(year), dat))
  
  m <- tapply(dat$SLOPE,as.numeric(dat$year), mean, na.rm = T)
  s <-  tapply(dat$SLOPE,as.numeric(dat$year), sd, na.rm = T)
  n <- table(as.numeric(dat$year))
  
  
  plot(msr ~m, type = "n")
  text(m, msr, labels = names(m))
  add.stats(f <- lm(msr ~ m))
  
}

#investigate role of SRali  effect  # always positive
par(mfrow = c(8,3), xpd = T, mar = c(3,3,2,1))
for (sp in impsp){
  dat <- databp[ databp$SpeciesCode == sp,]
  # barplot(table(dat$year,dat$abun) , beside = T )
  
  boxplot( SRali ~ as.numeric(year), dat)
  add.stats(f <- lm(abun ~ as.numeric(year), dat))
  
  msr <- tapply(dat$SRnat,as.numeric(dat$year), mean, na.rm = T)
  ssr <-  tapply(dat$SRnat,as.numeric(dat$year), sd, na.rm = T)
  
  plot(as.numeric(names(msr)), msr, ylim = c(-5,20), xlab = "months", ylab = "SRnat")
  segments(as.numeric(names(msr)), msr+ssr,as.numeric(names(msr)), msr -ssr)
  add.stats(f <- lm(SRnat ~ as.numeric(year), dat))
  
  m <- tapply(dat$SRali,as.numeric(dat$year), mean, na.rm = T)
  s <-  tapply(dat$SRali,as.numeric(dat$year), sd, na.rm = T)
  n <- table(as.numeric(dat$year))
  
  
  plot(msr ~m, type = "n")
  text(m, msr, labels = names(m))
  add.stats(f <- lm(msr ~ m))
  
}



### TEST GLM best regressions for all focal species:

glm.fits <- list()
sel.glm <- list()
for (i in 1:length(rownames(glmSRnat.overall[[1]]))){
 sp <- rownames(glmSRnat.overall[[1]])[i]
   dat <- databp[ databp$SpeciesCode == sp & databp$PlotName %in% unimprovedgrasslands,]
  f.bic <- bestglm( dat[c("abun","DOY.t","year","SLOPE","DEM_10", "Northern", "SRali", "SRnat")],family = poisson, IC = "BIC",RequireFullEnumerationQ = TRUE)
 print(f.bic$BestModels)
 glm.fits [[i]] <- f.bic$BestModel
 sel.glm  [[i]] <- f.bic$BestModels[1,]
 names(glm.fits)[i] <- sp
}

result <- do.call(rbind,sel.glm)
rownames(result) <- names(glm.fits)

impsp %in% rownames(result)[result$abun ==TRUE]  ## All significant species remain except ANTODO because non-monotonous

