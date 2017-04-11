### Model explaining gamma richness loss

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

# variation of gamma loss with alpha loss
# looking for a statistic capturing the slope of this trend, to be then related to some betadiversity parameter

# ES above abundance class
g = as.numeric(unlist(gamma.above.nat$gamma.above)) - as.numeric(unlist(gamma.above.nat$null.above))
a = as.numeric(unlist(alpha.above.trend$alpha.above)) - as.numeric(unlist(alpha.above.trend$null.above))
f <- lm(g ~ a)

# tmp = data.frame(a=a, g=g, sp = as.factor(rep(rownames(gamma.above.tre nd$gamma.above),6)))
# nlme(g ~ SSasymp(a , sp) , data=tmp, fixed = g ~ a,random =g ~ sp )

new = data.frame(a= seq(-10, 10, by=0.1))
pf <- predict(f,newdata=new,interval = "confidence")


plot(a,g, col="darkgrey", cex =0.6,
     ylim=c(-max(abs(g), na.rm=T), max(abs(g), na.rm=T)),
     xlim=c(-max(abs(a), na.rm=T), max(abs(a), na.rm=T)),
     xlab ="change in mean alpha richness",ylab ="change in gamma richness" )
abline(h=0, v=0, col="grey")
abline(f)
# lines(new$a,pf[,"lwr"], col="grey")
# lines(new$a,pf[,"upr"], col="grey")

col.impsp <-  c("grey", "black")[rownames(gamma.above.nat$gamma.above) %in% impsp +1]
names(col.impsp) <- rownames(gamma.trend.nat$obs)


# for (k in 1:length(rownames(gamma.above.nat$gamma.above))) {
for (k in 1:length(impsp)) {
  sp<- impsp[k]
  # sp<- rownames(gamma.above.nat$gamma.above) [k]
  g.sp <- as.numeric(gamma.above.nat$gamma.above[sp,] - gamma.above.nat$null.above[sp,])
  a.sp <- as.numeric(alpha.above.trend$alpha.above[sp,] - alpha.above.trend$null.above[sp,])

  col.sig.gam <-  c(NA, "black") [(gamma.above.nat$P.above[sp,]<0.05) +1]


  col.sig.gam <-  rep(NA, 6)
  col.sig.gam[glmSRnat.overall$impact.spread[sp,]$th.CI] <- "darkgrey"
  col.sig.gam[(gamma.above.nat$P.above[sp,]<0.05) &  col.sig.gam == "darkgrey"] <- "black"


  cex.sig.gam <-  rep(0.6, 6)
  cex.sig.gam[glmSRnat.overall$impact.spread[sp,]$th.CI] <- 1



  points(x = a.sp, y = g.sp, col = col.impsp[sp], bg= col.sig.gam, cex = cex.sig.gam,pch = 21)


  f<-  lm(g.sp ~ a.sp)
  pf <- data.frame(a = na.omit(a.sp), p = predict(f))
  pf[order(pf$a),]
  # lines(pf, col = col.impsp[sp])
  # lines(x = a.sp, y = g.sp, col = col.impsp[sp])

}


# SES above abundance class
g = as.numeric(unlist(gamma.above.nat$z.above))
a = as.numeric(unlist(alpha.above.trend$z.above))
f <- lm(g ~ a)

# tmp = data.frame(a=a, g=g, sp = as.factor(rep(rownames(gamma.above.nat$gamma.above),6)))
# nlme(g ~ SSasymp(a , sp) , data=tmp, fixed = g ~ a,random =g ~ sp )

new = data.frame(a= seq(-10, 10, by=0.1))
pf <- predict(f,newdata=new,interval = "confidence")

col.sp<-  c("grey", "red")[rep(rownames(gamma.above.nat$gamma.above),6) %in% aliens +1]
names(col.sp) <- rep(rownames(gamma.above.nat$gamma.above),6)



plot(a,g,cex =0.6, col = col.sp,
     ylim=c(-max(abs(g), na.rm=T), max(abs(g), na.rm=T)),
     xlim=c(-max(abs(a), na.rm=T), max(abs(a), na.rm=T)),
     xlab ="change in alpha richness",ylab ="change in gamma richness" )
abline(h=0, v=0, col="grey")
abline(f)
lines(new$a,pf[,"lwr"], col="grey")
lines(new$a,pf[,"upr"], col="grey")

col.impsp <-  c("grey", "black")[rownames(gamma.above.nat$gamma.above) %in% impsp +1]
names(col.impsp) <- rownames(gamma.trend.nat$obs)


# for (k in 1:length(rownames(gamma.above.nat$gamma.above))) {
for (k in 1:length(impsp)) {
  sp<- impsp[k]
  # sp<- rownames(gamma.above.nat$gamma.above) [k]
  g.sp <- as.numeric(gamma.above.nat$z.above[sp,] )
  a.sp <- as.numeric(alpha.above.trend$z.above[sp,])

  col.sig.gam <-  c(NA, "black") [(gamma.above.nat$P.above[sp,]<0.05) +1]


  col.sig.gam <-  rep(NA, 6)
  col.sig.gam[glmSRnat.overall$impact.spread[sp,]$th.CI] <- "black"

  cex.sig.gam <-  rep(NA, 6)
  cex.sig.gam[glmSRnat.overall$impact.spread[sp,]$th.CI] <- "black"


  points(x = a.sp, y = g.sp, col = col.impsp[sp], bg= col.sig.gam, pch = 21)


  f<-  lm(g.sp ~ a.sp)
  pf <- data.frame(a = na.omit(a.sp), p = predict(f))
  pf[order(pf$a),]
  # lines(pf, col = col.impsp[sp])
  # lines(x = a.sp, y = g.sp, col = col.impsp[sp])

}


# ES at abundance class
g = as.numeric(unlist(gamma.trend.nat$obs)) - as.numeric(unlist(gamma.trend.nat$mean.null))
a = as.numeric(unlist(alpha.trend.nat$obs)) - as.numeric(unlist(alpha.trend.nat$mean.null))
f <- lm(g ~ a)

new = data.frame(a=seq(-10, 10, by=0.1))
pf <- predict(f,newdata=new,interval = "prediction")

plot(a,g, type ="p", col= "grey", cex =0.6)
abline(h=0, v =0, col= "grey")
abline(f)
lines(new$a,pf[,"lwr"], lty="dotted")
lines(new$a,pf[,"upr"], lty="dotted")

col.impsp <-  c("grey", "black")[rownames(gamma.trend.nat$obs) %in% impsp +1]
names(col.impsp) <- rownames(gamma.trend.nat$obs)
K_slopes <- data.frame(matrix(NA, nrow = 79, ncol =4))
rownames(K_slopes) <- rownames(gamma.trend.nat$obs)
colnames(K_slopes) <- c("slope","r2","df","P")

# for (k in 1:length(rownames(gamma.trend.nat$obs))) {
for (k in 1:length(impsp)) {
  # sp<- rownames(gamma.trend.nat$obs) [k]
  sp<- impsp[k]

  g.sp <- as.numeric(gamma.trend.nat$obs[sp,] - gamma.trend.nat$mean.null[sp,])
  a.sp <- as.numeric(alpha.trend.nat$obs[sp,] - alpha.trend.nat$mean.null[sp,])

  col.sig.gam <-  rep(NA, 6)
  col.sig.gam[min(which((gamma.trend.nat$P.gamma[sp,]<0.05)), na.rm=T)] <- "black"

  col.sig.gam <-  rep(NA, 6)
  col.sig.gam[glmSRnat.overall$impact.spread[sp,]$th.CI] <- "black"

  points(x = a.sp, y = g.sp, col = col.impsp[sp], bg= col.sig.gam, pch = 21, cex = 0.6)

  f<-  lm(g.sp ~ a.sp)
  pf <- data.frame(a = na.omit(a.sp), p = predict(f))
  pf[order(pf$a),]
  # lines(pf, col = col.impsp[sp])
  # lines(x = a.sp, y = g.sp, col = col.impsp[k])
  sf = summary(f)
  K_slopes[k,] <- c(  slope <- f$coefficients[2], sf$adj.r.squared, sf$df[2], sf$coefficients[2,4] )
}


### correlation of slopes to a beta-diversity index :
beta.a = table.div.part[impsp,]$above.obs_sorensen
beta.b = table.div.part[impsp,]$below.obs_sorensen
alpha.a = table.div.part[impsp,]$aRc
alpha.b = table.div.part[impsp,]$aRnull

gamma.a = table.div.part[impsp,]$GRc
gamma.b = table.div.part[impsp,]$GRnull

est <- (beta.a - beta.b)/ (alpha.a - alpha.b) + beta.a  ## according to my fancy calculations...

K.th <- (gamma.a - gamma.b)/ (alpha.a - alpha.b)  ## direct calculation of ratio of loss : appears unrelated to the slope!

plot(K_slopes[impsp,]$slope,K.th)
plot(K_slopes[impsp,]$slope,est) # not at all related to the slope ! theory fails...
plot(est, K.th)

plot(K.th, table.div.part[impsp,]$below.obs_turnover)



############### gamma richness above threshold with NM2 = cumulated decline
par(mfrow = c(4,4), mar=c(0,0,1,1), oma=c(6,6,2,1))

list.data <-gamma.above.nat.betarare
gamma.loss <- as.data.frame(matrix(NA, nrow = length(impsp), ncol =6))
rownames(gamma.loss) <- impsp

ylim=c(-1,2)
ylim=c(0,250)
for (i in 1 : length(impsp)){

  if (i ==4) plot.new()
  sp <- impsp[i]
  # identify significant classes in black
  cols <- c(NA, "red") [ (list.data$P.above[sp,] <=0.025) +1]
  cols [ (list.data$P.above[sp,] >=0.975)] <- "black"

  n <- list.data$nplot.above[sp,]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y =  list.data$gamma.below[sp,]
  y [ (1:6)[ n < 5 | is.na(n)] ]=NA
  y[1]=NA

  gamma.loss[sp,] = y

  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-1.2, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)

  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  abline(h=0,lty="dotted")

  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)

  # Add species name
  mtext(3, text=paste(species[sp, "Genus"], "\n",species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=-2.5, las = 1)

  #   # Y axis label
  #   if ( i %in% c(1,4,7)) {
  #     mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
  #   }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("Cumulated loss in gamma richness"), adj=0.5, line=3, las = 0, outer=T)


############### Standardized gamma richness above threshold with NM2 = cumulated decline corrected for sample size
par(mfrow = c(4,4), mar=c(0,0,1,1), oma=c(6,6,2,1))

list.data <-gamma.above.nat
gamma.loss <- as.data.frame(matrix(NA, nrow = length(impsp), ncol =6))
rownames(gamma.loss) <- impsp

ylim=c(-5,5)
for (i in 1 : length(impsp)){

  if (i ==4) plot.new()
  sp <- impsp[i]
  # identify significant classes in black
  cols <- c(NA, "red") [ (list.data$P.above[sp,] <=0.025) +1]
  cols [ (list.data$P.above[sp,] >=0.975)] <- "black"

  n <- list.data$nplot.above[sp,]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y =  list.data$z.above[sp,]
  y [ (1:6)[ n < 5 | is.na(n)] ]=NA
  y[1]=NA

  gamma.loss[sp,] = y

  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-5.5, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)

  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  abline(h=0,lty="dotted")

  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)

  # Add species name
  mtext(3, text=paste(species[sp, "Genus"], "\n",species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=-2.5, las = 1)

  #   # Y axis label
  #   if ( i %in% c(1,4,7)) {
  #     mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
  #   }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("Cumulated change in gamma richness"), adj=0.5, line=3, las = 0, outer=T)


############### RAW DELTAgamma richness above threshold  + bootstrap + NM1 +NM0
par(mfrow = c(4,4), mar=c(0,0,1,1), oma=c(6,6,2,1))

list.data <-gamma.above.nat.betarare

ylim <- c(-250,250)
xlim <- c(-0.5, 5.5)
for (i in 1 : length(impsp)){

  if (i ==4) plot.new()
  sp <- impsp[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ ((list.data$P.delta[sp,] <=0.025)|(list.data$P.delta[sp,] >=0925)) +1]
  n <- list.data$nplot.above[sp,]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y =  list.data$deltagamma[sp,]
  y [ (1:6)[ n < 5 | is.na(n)] ]=NA
  # y[1]=NA


#   # null expectations
#   yn =  as.numeric(list.data$null.delta[sp,])
#   yn [ (1:6)[ n < 5 | is.na(n)] ]=NA
#   yn[1]=NA
#
#   sdn =  as.numeric(list.data$sdnull.delta[sp,])
#   sdn [ (1:6)[ n < 5 | is.na(n)] ]=NA
#   sdn[1]=NA


  # NULL CI NM0
  lownull0 = as.numeric(gamma.above.nat$lowCI.null.delta[sp,] )
  lownull0 [ (1:6)[ n < 5] ]=NA
  # lownull0[1]=NA
  hinull0 =  as.numeric(gamma.above.nat$hiCI.null.delta[sp,])
  hinull0 [ (1:6)[ n < 5] ]=NA
  # hinull0[1]=NA

  # NULL CI NM1
     lownull1 = as.numeric(list.data$lowCI.null.delta[sp,] )
     lownull1 [ (1:6)[ n < 5] ]=NA
     # lownull1[1]=NA
     hinull1 =  as.numeric(list.data$hiCI.null.delta[sp,])
     hinull1 [ (1:6)[ n < 5] ]=NA
     # hinull1[1]=NA


  # Bootstrapped CI
  lowCI = as.numeric(list.data$lowCI.boot.delta[sp,] )
  lowCI [ (1:6)[ n < 5 | is.na(n)] ]=NA
  # lowCI[1]=NA
  hiCI =  as.numeric(list.data$hiCI.boot.delta[sp,])
  hiCI [ (1:6)[ n < 5 | is.na(n)]]=NA
  # hiCI[1]=NA


  # plot background
  plot(x,y,ylim=ylim, xlim = xlim, las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-280, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)

  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  abline(h=0,lty="dotted")


  #plot NULL CI NM0
  cinull0 <-c(lownull0, hinull0[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull0)] <-NA
  polygon(na.omit(cix), na.omit(cinull0), bg="grey", col="grey")

  #plot NULL CI NM1
  cinull1 <-c(lownull1, hinull1[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull1)] <-NA
  polygon(na.omit(cix), na.omit(cinull1), bg="#FF000080", col="#FF000080")

  # plot bootstrap CI
  segments(x,lowCI,x, hiCI)

  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = xlim, xaxt = "n", yaxt="n", ann=F)


  # Add species name
  mtext(3, text=paste(species[sp, "Genus"], "\n",species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=-2.5, las = 1)

  #   # Y axis label
  #   if ( i %in% c(1,4,7)) {
  #     mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
  #   }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("Cumulated change in gamma richness"), adj=0.5, line=3, las = 0, outer=T)
mtext(3, text=c("NM1 sp resampling"), adj=0.5, line=0.51, las = 1, outer=T)

############### standardized effect size of DELTAgamma richness above threshold  + bootstrap + NM1
par(mfrow = c(3,4), mar=c(0,0,2,1), oma=c(6,6,4,1))

list.data <-gamma.above.nat

ylim=c(-50,50)
for (i in 1 : length(impsp)){

  if (i ==4) plot.new()
  sp <- impsp[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ (list.data$P.delta[sp,] <=0.025) +1]
  n <- list.data$nplot.above[sp,]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y =  list.data$deltagamma[sp,]  - list.data$null.delta[sp,]
  y [ (1:6)[ n < 5 | is.na(n)] ]=NA
  y[1]=NA

  # NULL CI NM1
  #    lowCI = as.numeric(list.data$lowCI.boot.delta[sp,] )
  #    lowCI [ (1:6)[ n < 5] ]=NA
  #    lowCI[1]=NA
  #    hiCI =  as.numeric(list.data$hiCI.boot.delta[sp,])
  #    hiCI [ (1:6)[ n < 5] ]=NA
  #    hiCI[1]=NA


  # Bootstrapped CI
  lowCI = as.numeric(list.data$lowCI.boot.delta[sp,] )
  lowCI [ (1:6)[ n < 5] ]=NA
  lowCI[1]=NA
  hiCI =  as.numeric(list.data$hiCI.boot.delta[sp,])
  hiCI [ (1:6)[ n < 5] ]=NA
  hiCI[1]=NA


  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-9, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)

  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  abline(h=0,lty="dotted")
  abline(h=0,lty="dotted")
  #
  #   #plot NULL CI NM1
  #   cinull1 <-c(lownull1, hinull1[6:1])
  #   cix <- c(x,x[6:1])
  #   cix[is.na(cinull1)] <-NA
  #   polygon(na.omit(cix), na.omit(cinull1), bg="#FF000080", col="#FF000080")

  # plot bootstrap CI
  # segments(x,lowCI,x, hiCI)

  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)

  # Add species name
  mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)

  #   # Y axis label
  #   if ( i %in% c(1,4,7)) {
  #     mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
  #   }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("SES (Difference in native gamma richness above Ac) "), adj=0.5, line=3, las = 0, outer=T)
mtext(3, text=c("NM1 sp resampling"), adj=0.5, line=0.51, las = 1, outer=T)



######### NM1 for gamma trends

par(mfrow = c(4,4), mar=c(0,0,2,1), oma=c(6,6,2,1))
ylim = c(-50,60)
xlim = c(0.5, 5.5)
sp.names = impsp

model <- gamma.trend.nat.beta$
for (i in 1 : length(sp.names)){
  sp <- sp.names[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ ((model$P.gamma[sp,] <=0.025)| (model$P.gamma[sp,] >=0.975)) +1]
  cols[1]=NA
  n <- as.numeric(glmSRnat.overall$n.obs[sp,])[1:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y = (model$obs[sp,] - model$mean.null[sp,])
  y [ (1:6)[ n < 5] ]=NA
  y[1]=NA

  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-57, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)

  axis(2, tcl= 0.1, label=F,  mgp=c(1,0.5,0), las=1)
  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  abline(h=0,lty="dotted")
  abline(v =glmSRnat.overall$impact.spread[sp,threshold]-1, lty="dashed")


  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)

  # Add species name
  mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)

  # Y axis label
  if ( i %in%  c(1,4,8)){
    mtext(2, text="ES",adj=0.5, cex=0.8, line=1.5, las = 0)
  }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4.5, las = 1, outer=T)
mtext(2, text=c("Native gamma richness(obs - null)"), adj=0.5, line=4, las = 0, outer=T)
mtext(3, text=c("NM1 sp resampling beta"), adj=0.5, line=0.51, las = 1, outer=T)



######### RAW VALUES gamma trends with different null model represented

par(mfrow = c(4,4), mar=c(0,0,2,1), oma=c(6,6,4,1))
ylim = c(0,200)
xlim = c(-0.5, 5.5)
sp.names = impsp
for (i in 1 : length(sp.names)){

  if (i ==4) {
    plot.new()
    legend( "topleft", legend = c("permute all", "permute with rare", "resample sp", "resample in rare"),
               fill = c("blue", "grey", "#FF000080","#EF9D5280"))
  }
  sp <- sp.names[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ ((gamma.trend.nat$P.gamma[sp,] <=0.025)| (gamma.trend.nat$P.gamma[sp,] >=0.975)) +1]
  cols[1]=NA
  n <- as.numeric(glmSRnat.overall$n.obs[sp,])[1:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y = gamma.trend.nat$obs[sp,]
  y [ (1:6)[ n < 5] ]=NA

  # NULL CI permute ALL
  lownull.all= as.numeric(gamma.trend.nat.permute.all$q025.null[sp,] )
  lownull.all [ (1:6)[ n < 5] ]=NA
  # lownull[1]=NA
  hinull.all =  as.numeric(gamma.trend.nat.permute.all$q975.null[sp,])
  hinull.all [ (1:6)[ n < 5] ]=NA
  # hinull[1]=NA


  # NULL CI NM0
  lownull= as.numeric(gamma.trend.nat$q025.null[sp,] )
  lownull [ (1:6)[ n < 5] ]=NA
  # lownull[1]=NA
  hinull =  as.numeric(gamma.trend.nat$q975.null[sp,])
  hinull [ (1:6)[ n < 5] ]=NA
  # hinull[1]=NA

  # NULL CI NM1
  lownull1= as.numeric(gamma.trend.nat.beta$q025.null[sp,] )
  lownull1 [ (1:6)[ n < 5] ]=NA
  # lownull1[1]=NA
  hinull1 =  as.numeric(gamma.trend.nat.beta$q975.null[sp,])
  hinull1 [ (1:6)[ n < 5] ]=NA
  # hinull1[1]=NA


  # NULL CI NM1
  lownull2= as.numeric(gamma.trend.nat.betarare$q025.null[sp,] )
  lownull2 [ (1:6)[ n < 5] ]=NA
  # lownull1[1]=NA
  hinull2 =  as.numeric(gamma.trend.nat.betarare$q975.null[sp,])
  hinull2 [ (1:6)[ n < 5] ]=NA
  # hinull1[1]=NA

  # plot background
  plot(x,y,ylim=ylim, xlim = xlim, las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-13, x = 0:5, labels= abclasses[1:6],  cex=1, srt=45, adj=1, xpd = NA)

  axis(2, tcl= 0.1, label=F,  mgp=c(1,0.5,0), las=1)
  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  abline(v =glmSRnat.overall$impact.spread[sp,threshold]-1, lty="dashed")

#  #plot NULL CI permute.all
#   cinull.all <-c(lownull.all, hinull.all[6:1])
#   cix <- c(x,x[6:1])
#   cix[is.na(cinull.all)] <-NA
#   polygon(na.omit(cix), na.omit(cinull.all), col="blue", border="blue")

  #plot NULL CI NM0
  cinull <-c(lownull, hinull[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull)] <-NA
  polygon(na.omit(cix), na.omit(cinull), bg="grey", col="grey")

  #plot NULL CI NM1
  cinull1 <-c(lownull1, hinull1[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull1)] <-NA
  polygon(na.omit(cix), na.omit(cinull1), bg="#FF0000", col="#FF0000")

#   #plot NULL CI NM1 resample in rare only
#   cinull2 <-c(lownull2, hinull2[6:1])
#   cix <- c(x,x[6:1])
#   cix[is.na(cinull2)] <-NA
#   polygon(na.omit(cix), na.omit(cinull2), col="#EF9D5280", border="#EF9D5280")


  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = xlim, xaxt = "n", yaxt="n", ann=F)

  # Add species name
  mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)

  # Y axis label
  if ( i %in%  c(1,4,8)){
    mtext(2, text="richness",adj=0.5, cex=0.8, line=1.7, las = 0)
  }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4.5, las = 1, outer=T)
mtext(2, text=c("Native gamma richness"), adj=0.5, line=4, las = 0, outer=T)


######### RAW VALUES DELTA gamma trends with both null model represented

par(mfrow = c(4,4), mar=c(0,0,2,1), oma=c(6,6,4,1))
ylim = c(-120,100)
xlim = c(-0.5, 5.5)
sp.names = impsp
for (i in 1 : length(sp.names)){


  if (i ==4) {
    plot.new()
    legend( "topleft", legend = c("permute all", "permute with rare", "resample sp", "resample in rare"),
            fill = c("blue", "grey", "#FF000080","#EF9D5280"))
  }

  sp <- sp.names[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ ((gamma.trend.nat.betarare$P.deltagamma[sp,] <=0.025)| (gamma.trend.nat.betarare$P.deltagamma[sp,] >=0.975)) +1]
  cols[1]=NA
  n <- as.numeric(glmSRnat.overall$n.obs[sp,])[1:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y = gamma.trend.nat$obs.deltagamma[sp,]
  y [ (1:6)[ n < 5] ]=NA

  # NULL CI permute ALL
  lownull.all= as.numeric(gamma.trend.nat.permute.all$q025.null.deltagamma[sp,] )
  lownull.all [ (1:6)[ n < 5] ]=NA
  # lownull[1]=NA
  hinull.all =  as.numeric(gamma.trend.nat.permute.all$q975.null.deltagamma[sp,])
  hinull.all [ (1:6)[ n < 5] ]=NA
  # hinull[1]=NA

  # NULL CI NM0
  lownull= as.numeric(gamma.trend.nat$q025.null.deltagamma[sp,] )
  lownull [ (1:6)[ n < 5] ]=NA
  # lownull[1]=NA
  hinull =  as.numeric(gamma.trend.nat$q975.null.deltagamma[sp,])
  hinull [ (1:6)[ n < 5] ]=NA
  # hinull[1]=NA

  # NULL CI NM1
  lownull1= as.numeric(gamma.trend.nat.beta$q025.null.deltagamma[sp,] )
  lownull1 [ (1:6)[ n < 5] ]=NA
  # lownull1[1]=NA
  hinull1 =  as.numeric(gamma.trend.nat.beta$q975.null.deltagamma[sp,])
  hinull1 [ (1:6)[ n < 5] ]=NA
  # hinull1[1]=NA

  # NULL CI NM1 reample in rare only
  lownull2= as.numeric(gamma.trend.nat.betarare$q025.null.deltagamma[sp,] )
  lownull2 [ (1:6)[ n < 5] ]=NA
  # lownull1[1]=NA
  hinull2 =  as.numeric(gamma.trend.nat.betarare$q975.null.deltagamma[sp,])
  hinull2 [ (1:6)[ n < 5] ]=NA
  # hinull1[1]=NA

  # plot background
  plot(x,y,ylim=ylim, xlim = xlim, las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-130, x = 0:5, labels= abclasses[1:6],  cex=1, srt=45, adj=1, xpd = NA)

  axis(2, tcl= 0.1, label=F,  mgp=c(1,0.5,0), las=1)
  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  abline(v =glmSRnat.overall$impact.spread[sp,threshold]-1, lty="dashed")

  #plot NULL CI permute.all
  cinull.all <-c(lownull.all, hinull.all[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull.all)] <-NA
  polygon(na.omit(cix), na.omit(cinull.all), col="blue", border="blue")


  #plot NULL CI NM0
  cinull <-c(lownull, hinull[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull)] <-NA
  polygon(na.omit(cix), na.omit(cinull), col="grey", border="grey")

  #plot NULL CI NM1
  cinull1 <-c(lownull1, hinull1[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull1)] <-NA
  polygon(na.omit(cix), na.omit(cinull1), col="#FF000080", border="#FF000080")


  #plot NULL CI NM1 resample in rare only
  cinull2 <-c(lownull2, hinull2[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull2)] <-NA
  polygon(na.omit(cix), na.omit(cinull2), col="#EF9D5280", border="#EF9D5280")

  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = xlim, xaxt = "n", yaxt="n", ann=F)

  # Add species name
  mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)

  # Y axis label
  if ( i %in%  c(1,4,8)){
    mtext(2, text="richness",adj=0.5, cex=0.8, line=1.7, las = 0)
  }
}
mtext(3, text=c("grey : NM0 permute.rare\nred : NM1 sp resampling.all"), adj=0.5, line=0.51, las = 1, outer=T)

mtext(1, text=c("Abundance class"), adj=0.5, line=4.5, las = 1, outer=T)
mtext(2, text=c("Difference in Native gamma richness"), adj=0.5, line=4, las = 0, outer=T)

