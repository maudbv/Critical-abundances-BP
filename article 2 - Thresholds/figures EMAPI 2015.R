# Figure Article thresholds of impact

### import data
# load("saved Rdata/article 2 - threshold/article threshold 1.2.Rdata")

abclasses= c("Rare" ,"Occasional",  "Frequent", "Common", "Abundant","Dominant")
db <- databp[databp$PlotName %in% realgrasslands,]

#########  Details per significant species ##############
x11()

threshold = "th.CI"
sp<-"LOLPER"
Z= db[db$SpeciesCode==sp,]
M <- glmSRnat.overall$impact.spread
variable <- "SRnat"
es <- as.numeric(M$mean.values[sp,])

#symbol colors
col = rep("grey",length(Z$abun))
# col[ Z$abun == as.numeric( M$thresh[sp,threshold])]="black"  # identify the threshold

# jitter points
Z$abun.jit = jitter(Z$abun, amount = 0.4)

#plotting
plot(as.formula(paste(variable, " ~ abun.jit")), data=Z, xlim=c(0.5,6.5),
     type="n",xaxt = "n", yaxt="n", ann=F, bty = "l")

points(as.formula(paste(variable, " ~ abun.jit")), data=Z, xlim=c(0.5,6.5),
       col="grey", bg=col,pch=21,  xaxt = "n", yaxt="n", ann=F, cex = 0.8)

axis(1,at = 1:6, label = abclasses, tcl= 0.1,mgp=c(1,0.5,0),las=2)
axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)


points(1:6,es, xlim=c(0.5,6.5),col="black", bg="black" ,pch=21, cex = 1)

abline(h = mean(Z[Z$abun == 1, variable], na.rm=T), lty="dotted")
segments(1:6,es[1],1:6,es, xlim=c(0.5,6.5),col="black")


## with GLM effect size

es <- as.numeric(M$est[sp,])
n <- as.numeric(M$n.obs[sp,])[2:6]

low <- as.numeric(M$CIlow[sp,])
hi <- as.numeric(M$CIhi[sp,])

## correct for small sample sizes
small <- es
small[which(n>=3)] <- NA
es[which(n<3)] <- NA
low[n<3] <- NA
hi[n<3] <- NA

## ylimns
if ( !all(is.na(c(hi, low)))) lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
if ( all(is.na(c(hi, low)))) lims <-  max(c(abs(es)), na.rm=T) +0.02
ylim=c(-lims, lims)

# symbol color
M$thresh <-M$impact.spread
col <- rep("white",6)
col[ as.numeric( M$thresh[sp,threshold])] <- "black"  # identify the threshold

plot(1:5,  rep(0,5), ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F)
axis(1,at = 1:5, label = abclasses[2:6], tcl= 0.1,mgp=c(1,0.5,0),las=2)
axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
abline(h=0,lty="dotted")

## Add bootstrapped CI
if ( !all(is.na(c(hi, low)))) {
  arrows(1:5,low,   1:5,hi, lwd=1, code =3, length=0.05, angle=90)
}
points(1:5,es, bg = col[2:6], pch=21)
mtext(side=3, line=-1.2, at=1:5,text=n, cex=0.7, las=1)


############## Barplot of frequencies for Alien richness only #############

n <- length(unique(glmSRnat.sum$class.summary$group))
par(mfcol=c(1,1), oma=c(1,1,1,1), mar=c(4,1,1,1), las=1)
sum.df = glmSRnat.sum$class.summary
  S <- as.data.frame(sum.df[sum.df$group == unique(sum.df$group)[2],])


  barplot(S$nb.sp, ylim=c(0,max(40, S$nb.sp)),
          col= "grey",  border= NA, axes=F)

  par(new=T)
  barplot(S$freq.negative.sig, ylim=c(0,max(40, S$nb.sp)),
          col= "grey50",  border= NA, axes=F)

  par(new=T)
  b <- barplot(S$freq.thr, col="orangered" , ylim=c(0,max(40, S$nb.sp)), border= NA, axes=F)
  axis(2, tcl = -0.3, cex.axis = 0.8, mgp = c(1,0.5,0), las = 1)
  text(y=-3, x = b, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)



############## Proportion of critical abundances for Aliens #############
x11()
n <- length(unique(glmSRnat.sum$class.summary$group))
par(mfcol=c(1,2), oma=c(5,6,4,1), mar=c(1,1,1,2), las=1)
sum.df = glmSRnat.sum$class.summary
 i=2
  S <- as.data.frame(sum.df[sum.df$group == unique(sum.df$group)[i],])


  barplot(rep(1,5), ylim=c(0,1),
          col= "grey",  border= NA, axes=F)

  par(new=T)
  barplot(S$prop.negative.sig, ylim=c(0,1),
          col= "grey50",  border= NA, axes=F)

  par(new=T)
  b <- barplot(S$prop.thr, col="orangered" , ylim=c(0,1), border= NA, axes=F)
  axis(2, tcl = -0.3, cex.axis = 0.8, mgp = c(1,0.5,0), las = 1)
  text(y=-0.05, x = b, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)



x11()
n <- length(unique(glmSRali.sum$class.summary$group))
par(mfcol=c(1,2), oma=c(5,6,4,1), mar=c(1,1,1,2), las=1)
sum.df = glmSRali.sum$class.summary
for (i in 1:n) {
  S <- as.data.frame(sum.df[sum.df$group == unique(sum.df$group)[i],])


  barplot(rep(1,5), ylim=c(0,1),
          col= "palegreen3",  border= NA, axes=F)

  par(new=T)
  barplot(S$prop.negative, ylim=c(0,1),
          col= "grey50",  border= NA, axes=F)

  par(new=T)
  b <- barplot(S$prop.thr, col="orangered" , ylim=c(0,1), border= NA, axes=F)
  axis(2, tcl = -0.3, cex.axis = 0.8, mgp = c(1,0.5,0), las = 1)
  text(y=-0.05, x = b, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)
}

##########  ALPHA trends effect sizes for ALL IMPSP ################
threshold ="th.CI"
sel <- impsp


par(mfrow = c(3,4), mar=c(0,0,2,2), oma=c(6,6,2,3))

M <- glmSRnat.overall
ylim=c(-3,3)
### Loop on selected species
for (i in 1:length(sel))  {
  if (i ==4) plot.new()
  sp <- sel[i]
  es <- as.numeric(M$est[sp,])
  n <- as.numeric(M$n.obs[sp,])[2:6]

  low <- as.numeric(M$CIlow[sp,])
  hi <- as.numeric(M$CIhi[sp,])
  M$thresh <-M$impact.spread

  ## correct for small sample sizes
  small <- es
  small[which(n>=3)] <- NA
  es[which(n<3)] <- NA
  low[n<3] <- NA
  hi[n<3] <- NA

  if ( !all(is.na(c(hi, low)))) lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
  if ( all(is.na(c(hi, low)))) lims <-  max(c(abs(es)), na.rm=T) +0.02

  col <- rep("white",6)
  col[(1:6)>= as.numeric( M$thresh[sp,threshold])] <- "black"  # above threshold
  col[(2:6)[hi>0]] <- NA  # robust negative coef
  col[(2:6)[n<5]] <- NA  # sufficient data points


  plot(1:5,  rep(0,5),ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8:11)) text(y=-3.6, x = 1:5, labels= abclasses[2:6],  cex=0.9, srt=45, adj=1, xpd = NA)

  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  abline(h=0,lty="dotted")
  abline(v = M$thresh[sp,threshold]-1, lty="dashed")

  if (!all(is.na(small))) points(1:5, small, pch=".")

  ## Add bootstrapped CI
  if ( !all(is.na(c(hi, low)))) {
    arrows(1:5,low,   1:5,hi, lwd=1, code =3, length=0.05, angle=90)
  }

  # draw the points and lines
  par(new=T)
  plot(1:5,es, bg = col[2:6], pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)

  # Add species name
  mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)


  # Y axis label
  if ( i %in% c(1,4,8)) {
    mtext(2, text="GLM coefficient", adj=0.5, cex=0.7, line=1.5, las = 0)
  }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("Effect on native alpha richness"), adj=0.5, line=4, las = 0, outer=T)

#########  Alpha richness ES for example species ##############
x11()
quartz()
threshold = "th.CI"
sp<-"ACHMIL"
sp<-"DACGLO"
sp<-"LOLPER"
Z= db[db$SpeciesCode==sp,]
M <- glmSRnat.overall
variable <- "SRnat"
es <- as.numeric(M$mean.values[sp,])

## with GLM effect size

es <- as.numeric(M$est[sp,])
n <- as.numeric(M$n.obs[sp,])[2:6]

low <- as.numeric(M$CIlow[sp,])
hi <- as.numeric(M$CIhi[sp,])

## correct for small sample sizes
small <- es
small[which(n>=3)] <- NA
es[which(n<3)] <- NA
low[n<3] <- NA
hi[n<3] <- NA

## ylimns
if ( !all(is.na(c(hi, low)))) lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
if ( all(is.na(c(hi, low)))) lims <-  max(c(abs(es)), na.rm=T) +0.02
ylim=c(-lims, lims)

# symbol color
col <- rep("white",6)
col[(1:6)>= as.numeric( M$impact.spread[sp,threshold])] <- "black"  # above threshold
col[(2:6)[hi>0]] <- NA  # robust negative coef
col[(2:6)[n<5]] <- NA  # sufficient data points



plot(1:5,  rep(0,5), ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F)
axis(1,at = 1:5, label = abclasses[2:6], tcl= 0.1,mgp=c(1,0.5,0),las=2)
axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
abline(h=0,lty="dotted")

## Add bootstrapped CI
if ( !all(is.na(c(hi, low)))) {
  arrows(1:5,low,   1:5,hi, lwd=1, code =3, length=0.05, angle=90)
}
points(1:5,es, bg = col[2:6], pch=21)



####### plotting impact spread ###########

db=databp[databp$PlotName %in% realgrasslands,]
targets <- c("ACHMIL", "DACGLO")

par(mfrow=c(1,3), mar=c(1,1,1,1))



for( i in targets) {
  plot.alien=as.character(db[which(db$abun%in%c(1,2,3,4,5,6) & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])
  th= glmSRnat.overall$impact.spread[i, "th.CI"]
  es <- as.numeric(glmSRnat.overall$est[i,])
  n <- as.numeric(glmSRnat.overall$n.obs[i,])[2:6]
  es[which(n<5)] <- NA
  es = which(!is.na(es))+1
  m <- max(es, na.rm=T)
  max.impact=as.character(db[which(db$abun>= m & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])
  th.impact=as.character(db[which(db$abun>= th & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])

  plot(study_area, col="white", border = "grey60")
  plot(envplot[plot.alien, ],pch=22, cex = 0.9, col="tan",bg="tan", add=T)
  plot(envplot[th.impact, ],pch=22, cex = 0.9, col="sienna2", bg="sienna2", add=T)
  plot(envplot[max.impact, ],pch=22,  cex = 0.9, col="sienna4", bg="sienna4", add=T)


  print(length(plot.alien))
  print(length(max.impact))
  print( length(th.impact))
  mtext(3, text =sub("_", " ",species[i, "tip"]), line= -1, font=3, cex = 0.7)
}

######### RAW VALUES gamma trends with different null model represented #####

par(mfrow = c(3,4), mar=c(0,0,2,1), oma=c(6,6,4,1))
ylim = c(0,200)
xlim = c(-0.5, 5.5)

  sp <- "LOLPER"
  # identify significant classes in black
  cols <- c(NA, "black") [ ((gamma.trend.nat.permute.all$P.gamma[sp,] <=0.025)| (gamma.trend.nat.permute.all$P.gamma[sp,] >=0.975)) +1]
  cols[1]=NA
  n <- as.numeric(glmSRnat.overall$n.obs[sp,])[1:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y = gamma.trend.nat.betarare$obs[sp,]
  y [ (1:6)[ n < 5] ]=NA


  # NULL CI permute ALL
  lownull.all= as.numeric(gamma.trend.nat.permute.all$q025.null[sp,] )
  lownull.all [ (1:6)[ n < 5] ]=NA
  # lownull[1]=NA
  hinull.all =  as.numeric(gamma.trend.nat.permute.all$q975.null[sp,])
  hinull.all [ (1:6)[ n < 5] ]=NA
  # hinull[1]=NA

  # NULL CI NM1
  lownull1= as.numeric(gamma.trend.nat.beta$q025.null[sp,] )
  lownull1 [ (1:6)[ n < 5] ]=NA
  # lownull1[1]=NA
  hinull1 =  as.numeric(gamma.trend.nat.beta$q975.null[sp,])
  hinull1 [ (1:6)[ n < 5] ]=NA
  # hinull1[1]=NA


  #
  #   # NULL CI NM1
  #   lownull2= as.numeric(gamma.trend.nat.betarare$q025.null[sp,] )
  #   lownull2 [ (1:6)[ n < 5] ]=NA
  #   # lownull1[1]=NA
  #   hinull2 =  as.numeric(gamma.trend.nat.betarare$q975.null[sp,])
  #   hinull2 [ (1:6)[ n < 5] ]=NA
  #   # hinull1[1]=NA

  # plot background
  plot(x,y,ylim=ylim, xlim = xlim, las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-13, x = 0:5, labels= abclasses[1:6],  cex=1, srt=45, adj=1, xpd = NA)

  axis(2, tcl= 0.1, label=F,  mgp=c(1,0.5,0), las=1)
  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  # abline(v =glmSRnat.overall$impact.spread[sp,threshold]-1, lty="dashed")

  #plot NULL CI permute.all
  cinull.all <-c(lownull.all, hinull.all[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull.all)] <-NA
  polygon(na.omit(cix), na.omit(cinull.all), col="grey", border="grey")

  #plot NULL CI NM1
  cinull1 <-c(lownull1, hinull1[6:1])
  cix <- c(x,x[6:1])
  cix[is.na(cinull1)] <-NA
  polygon(na.omit(cix), na.omit(cinull1), col="#FF000080", border="#FF000080")
  #
  #   #plot NULL CI NM1 resample in rare only
  #   cinull2 <-c(lownull2, hinull2[6:1])
  #   cix <- c(x,x[6:1])
  #   cix[is.na(cinull2)] <-NA
  #   polygon(na.omit(cix), na.omit(cinull2), col="#FF000080", border="#FF000080")


  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = xlim, xaxt = "n", yaxt="n", ann=F)

  # Add species name
  mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)

}

mtext(1, text=c("Abundance class"), adj=0.5, line=4.5, las = 1, outer=T)
mtext(2, text=c("Native gamma richness"), adj=0.5, line=4, las = 0, outer=T)

##### ES gamma trends AGAINST PERMUTE ALL ############

par(mfrow = c(3,4), mar=c(0,0,2,1), oma=c(6,6,2,1))
ylim = c(-50,35)
xlim = c(0.5, 5.5)

for (i in 1 : length(impsp)){

  sp <- impsp[i]
  if (i ==4) plot.new()
  # identify significant classes in black
  cols <- c(NA, "black") [ (gamma.trend.nat.permute.all$P.gamma[sp,] <=0.025 | gamma.trend.nat.permute.all$P.gamma[sp,] >=0.975) +1]
  n <- as.numeric(glmSRnat.overall$n.obs[sp,])[2:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y = (gamma.trend.nat.permute.all$obs[sp,] - gamma.trend.nat.permute.all$mean.null[sp,])
  # /gamma.trend.nat$sd.null[sp,]
  y [ (2:6)[ n < 5] ]=NA
  y[1]=NA


  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-56, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)

  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  abline(h=0,lty="dotted")
  abline(v =glmSRnat.overall$impact.spread[sp,threshold]-1, lty="dashed")

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
mtext(2, text=c("Effect size of native gamma richness"), adj=0.5, line=3, las = 0, outer=T)

##### ES gamma trends AGAINST PERMUTE ALL ############

par(mfrow = c(8,10), mar=c(0,0,2,1), oma=c(2,2,2,2))
ylim = c(-50,35)
xlim = c(0.5, 5.5)

for (i in 1 : length(rownames(glmSRnat.overall$n.obs))){

  sp <-rownames(glmSRnat.overall$n.obs)[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ (gamma.trend.nat.permute.all$P.gamma[sp,] <=0.025 | gamma.trend.nat.permute.all$P.gamma[sp,] >=0.975) +1]
  n <- as.numeric(glmSRnat.overall$n.obs[sp,])[2:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)

  # create y axis with the gamma richness Standardized effect size
  y = (gamma.trend.nat.permute.all$obs[sp,] - gamma.trend.nat.permute.all$mean.null[sp,])
  # /gamma.trend.nat$sd.null[sp,]
  y [ (2:6)[ n < 5] ]=NA
  y[1]=NA


  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  abline(h=0,lty="dotted")
  abline(v =glmSRnat.overall$impact.spread[sp,threshold]-1, lty="dashed")

  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)

  # Add species name
  mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.6, line=0.2, las = 1)

  #   # Y axis label
  #   if ( i %in% c(1,4,7)) {
  #     mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
  #   }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("Effect size of native gamma richness"), adj=0.5, line=3, las = 0, outer=T)


###  cORRELATIONS ##########

# raw values of loss in alpha vs. loss in gamma
x11()
dg<- table.div.part$deltagamma
da<-  table.div.part $aRc -table.div.part$aRo
cgam <- c("grey50", "black")[ (table.div.part$GRP.permute.all<=0.05) + 1]
plot(da, dg ,pch = 20, col=cgam , ann=F, type ="n")
text(x = da,y = dg,label = rownames(table.div.part) , col=cgam, cex =0.85)
f <- cor.test(da ,dg)
mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8)
mtext(1,text = expression(Delta*" alpha richness"), line=2.5, cex =0.85 )
mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)


### Pissible delta gamma above critical abudnances
dg<- table.div.part$deltagamma #raw
dg<- table.div.part$deltagamma -  table.div.part$delta.null.permute.all
dg<- table.div.part$GRc -  table.div.part$GRnull
# dg<- table.div.part$delta.z.permute.all
# dg<- table.div.part$delta.z.beta

cgam <- c("grey50", "black")[ (table.div.part$GRP.permute.all<=0.05) + 1]
# cgam <- c("grey50", "black")[ (table.div.part$GRP.beta<=0.05) + 1]


x11()
da<-  table.div.part $aRc -table.div.part$aRo
plot(da, dg ,pch = 20, col=cgam , ann=F, type ="n", las=1, xlim= c(-4.5, -1))
text(x = da,y = dg,label = rownames(table.div.part) , col=cgam, cex =0.85)
points(x = da,y = dg,pch = 20, col=cgam, cex =1)
f <- cor.test(da ,dg)
# mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8)
# mtext(1,text = expression(Delta*" alpha richness"), line=2.5, cex =0.85 )
# mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)

## threshold
thresh <-  glmSRnat.overall$impact.spread[rownames(table.div.part),  "th.CI"]
plot (thresh ,dg , pch = 20, col=cgam , ann=F, xlim= c(1.5,6.5), type = "n", las=1)
f <- cor.test(thresh ,dg)
# points(x = thresh,y = dg,pch = 20, col=cgam, cex =1)
text(x = thresh , y = dg,label =  rownames(table.div.part) ,col=cgam, cex =0.85 )
# mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8)
# mtext(1,text = "Critical abundance", line=2.5, cex =0.85 )
# mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)

### spread
spread <-  glmSRnat.overall$impact.spread[rownames(table.div.part),  "n.plot.impact"]
plot (spread,dg , pch = 20, col=cgam , ann=F, log = "x", xlim =c(2.5, 700), type= "n", las=1)
f <- cor.test(spread ,dg)
text(x = spread,y = dg,label = rownames(table.div.part) ,col=cgam, cex =0.85)
# mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8)
# mtext(1,text = "Number of plots > critical abundance", line=2.5, cex =0.85 )
# mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)

x11()
# prop of impact
spread <-  glmSRnat.overall$impact.spread[rownames(table.div.part),  "prop.plot.impact"]
plot (spread,dg , pch = 20, col=cgam , ann=F, xlim =c(0, 1), type="n")
f <- cor.test(spread ,dg)
text(x = spread, y = dg,label = rownames(table.div.part) ,col=cgam, cex =0.7 )
mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 0, cex=0.8)
mtext(1,text = "Prop. of plots > critical abundance", line=2.5, cex =0.85 )
mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)


###  alpha loss vs threshold
da<-  -( table.div.part $aRc - table.div.part$aRo )
thresh <-  glmSRnat.overall$impact.spread[rownames(table.div.part),  "th.CI"]
plot (thresh,da, pch = 20, col=cgam , ann=F,  type="n")
f <- cor.test(da,thresh, method = "spearman")
text(x = thresh,y = da, label = impsp ,col="black", cex =0.7 )
mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 0, cex=0.8)
mtext(1,text = "critical abundance", line=2.5, cex =0.85 )
mtext(2,text ="alpha loss", line=2.5  , las=0, cex =0.85)

