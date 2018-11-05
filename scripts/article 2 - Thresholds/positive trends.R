# Results for positive trends

# 2 alien species with positive trends:
 impsp.pos <- rownames(glmSRnat.overall$impact.spread[which(!is.na(glmSRnat.overall$impact.spread$pth.CI) & (rownames(glmSRnat.overall$impact.spread) %in% aliens)),])
# 

## With native focal species included:
# impsp.pos <- rownames(glmSRnat.overall$impact.spread[which(!is.na(glmSRnat.overall$impact.spread$pth.CI)),])

#### Extract result table for loss in native diversity above thresholds:   ######
sp.names = impsp.pos
out <- glmSRnat.overall$impact.spread
out <- out[sp.names,c("pth.CI", "prevalence", "n.plot.impact", "n.plot.dominant")]
out<- cbind(species = species[sp.names, "tip"], out)

out$pth.CI.SRali <- glmSRali.overall$impact.spread[sp.names,"pth.CI"]
out$th <- out$pth.CI
out$th[is.na(out$th)] <- out$pth.CI.SRali[is.na(out$th)]

out$aRc.old <- sapply(sp.names, FUN = function(m) glmSRnat.overall$mean[m,out[m, "th"]] )
out$aRo.old <- glmSRnat.overall$mean[sp.names, "C1"]

out$aRo <- sapply(sp.names, FUN = function(m) alpha.above.trend$alpha.below[m,out[m, "th"]] )
out$aRc <- sapply(sp.names, FUN = function(m) alpha.above.trend$alpha.above[m,out[m, "th"]] )
out$aRo.sd <- sapply(sp.names, FUN = function(m) alpha.above.trend$alpha.below.sd[m,out[m, "th"]] )
out$aRc.sd <- sapply(sp.names, FUN = function(m) alpha.above.trend$alpha.above.sd[m,out[m, "th"]] )

out$aRnull <- sapply(sp.names, FUN = function(m) alpha.above.trend$null.above[m,out[m, "th"]]  )
out$aRsdnull <- sapply(sp.names, FUN = function(m)  alpha.above.trend$sdnull.above[m,out[m, "th"]]  )
out$alpha.loss <- out$aRnull  - out$aRc
out$aR.P <- sapply(sp.names, FUN = function(m)  alpha.above.trend$P.above[m,out[m, "th"]] )

out$deltaalpha <- sapply(sp.names, FUN = function( m) alpha.above.trend$deltaalpha[m,out[m, "th"]]  )
out$deltaalpha.null <- sapply(sp.names, FUN = function(m)  alpha.above.trend$null.delta[m,out[m, "th"]]  )
out$deltaalpha.z.permute.all <- sapply(sp.names, FUN = function(m) {
  (alpha.above.trend$deltaalpha[m,out[m, "th"]] - alpha.above.trend$null.delta[m,out[m, "th"]])/alpha.above.trend$sdnull.delta[m,out[m, "th"]]
})

out$aR.t <- sapply(sp.names, FUN = function(m)  alpha.above.trend$t[m,out[m, "th"]] )
out$aR.df<- sapply(sp.names, FUN = function(m)  alpha.above.trend$df[m,out[m, "th"]] )
out$aR.Pt <- sapply(sp.names, FUN = function(m)  alpha.above.trend$P.t[m,out[m, "th"]] )

out$GRo <- sapply(sp.names, FUN = function(m) gamma.above.nat.permute.all$gamma.below[m,out[m, "th"]] )
out$GRc <- sapply(sp.names, FUN = function(m) gamma.above.nat.permute.all$gamma.above[m,out[m, "th"]] )
out$GRnull <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$null.above[m,out[m, "th"]]  )
out$GRsdnull <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$sdnull.above[m,out[m, "th"]]  )
out$gamma.loss <- out$GRnull - out$GRc

out$deltagamma <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$deltagamma[m,out[m, "th"]]  )
out$delta.z.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$z.delta[m,out[m, "th"]]  )
out$delta.null.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$null.delta[m,out[m, "th"]]  )
out$GRP.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$P.delta[m,out[m, "th"]] )

# out$deltagamma.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$deltagamma[m,out[m, "th"]]  )
# out$delta.z.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$z.delta[m,out[m, "th"]]  )
# out$delta.null.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$null.delta[m,out[m, "th"]]  )
# out$GRP.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$P.delta[m,out[m, "th"]] )

out$BRo <- sapply(sp.names , FUN = function(m) gamma.above.nat.permute.all$betap.below[m,out[m, "th"]] )
out$BRc <-sapply(sp.names, FUN = function(m) gamma.above.nat.permute.all$betap.above[m,out[m, "th"]] )

# out$z.beta.diff <-table.turnover[sp.names,"z.diff"]
# out$P.beta.diff <-table.turnover[sp.names,"p.diff"]

# final formatting of table
table.pos <- out[impsp.pos,]


table2.pos <- data.frame(Acrit = table.pos$pth.CI,
                     Presence = table.pos$prevalence,
                     above.Acrit = table.pos$n.plot.impact,
                     dominance = table.pos$n.plot.dominant,
                     alpha.below = round(table.pos$aRo, 1),
                     alpha.below.sd = round(table.pos$aRo.sd,1),
                     alpha.above = round(table.pos$aRc,1),
                     alpha.above.sd = round(table.pos$aRc.sd,1),
                     delta.alpha = table.pos$aRc -  table.pos$aRo,
                     perc.delta.alpha = round(((table.pos$aRc -  table.pos$aRo)/table.pos$aRo)*100, 1),
                     gamma.below = table.pos$GRo,
                     gamma.above = table.pos$GRc,
                     delta.gamma.c = table.pos$deltagamma - table.pos$delta.null.permute.all,
                     delta.gamma.P = 1 - round(table.pos$GRP.permute.all, 4)
)
rownames(table2.pos) <- table.pos$species


write.csv(table2.pos, file = "table2.positives.csv")

plot(table2.pos$Acrit, table2.pos$delta.alpha)






#### __________________________________________DEFINE ELEMENTS__________________________________________ ####
abclasses= c("Rare" ,"Occasional",  "Frequent", "Common", "Abundant","Dominant")
threshold ="pth.CI"
sel <- impsp.pos


#spatial distributions:
spread <- table.pos$n.plot.impact
prevalence <- table.pos$prevalence
dominance <- table.pos$n.plot.dominant
prop.impact <- spread/prevalence

# critical abundance level
critical.abun <- table.pos$th

# library(FactoMineR)
# tmp <- PCA(data.frame(cbind(dg0, dg1,dg, da,da1, dbeta), row.names= impsp ) )
# plot(tmp, choix = "ind")
# plot(tmp, choix = "var")
# plot(tmp, choix = "var", axes = c(1,3))

### Looking at GLM results across species and covariables
names(glmSRnat.overall)[6] <- "covar.tab"



##colour code For figures
# cgam <- c("grey50", "black")[ (table.pos$GRP.permute.all<=0.05) + 1]
cgam <- "black"


#### Figure S1: barplot of frequencies for positive trends:    #########

# extract frequency table:
effects <- glmSRnat.sum$class.summary
n <- length(unique(effects$group))

#draw barplot
par(mfrow = c(1,2), mar=c(5,3,7,3), las=1)
S <- effects[effects$group == "ALIEN:1",]
barplot(S$nb.sp, ylim=c(0,max(40, S$nb.sp)),col= "grey80",  border= NA, axes=F)
par(new=T)
barplot(S$freq.positive.above, ylim=c(0,max(40, S$nb.sp)),col= "grey60",  border= NA, axes=F)
par(new=T)
b <- barplot(S$freq.pthr, col="black" , ylim=c(0,max(40, S$nb.sp)), border= NA, axes=F)
axis(2, tcl = -0.3, cex.axis = 0.8, mgp = c(1,0.5,0), las = 1)
text(y=-1.5, x = b+0.3, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)
legend(x=0, y=60, bty="n", bg="white", 
       legend=c( "Total number of focal species",
                 "Number of positive effects",
                 "Number of positive critical abundances"),
       fill=c( "grey80","grey60","black"), border= c("grey90","grey60","black"), cex=0.8, xpd=NA,y.intersp =1)
mtext(text=c("a) Alien focal species"), side=3, outer=F, line=5)
mtext(text="number of species", side=2, outer=F, line=2, las=0, cex=1)
mtext(text="Abundance class", side=1, outer=F, line=3.5)

      


#draw barplot for natives
S <- effects[effects$group == "ALIEN:0",]
barplot(S$nb.sp, ylim=c(0,max(40, S$nb.sp)),col= "grey80",  border= NA, axes=F)
par(new=T)
barplot(S$freq.positive.above, ylim=c(0,max(40, S$nb.sp)),col= "grey60",  border= NA, axes=F)
par(new=T)
b <- barplot(S$freq.pthr, col="black" , ylim=c(0,max(40, S$nb.sp)), border= NA, axes=F)
axis(2, tcl = -0.3, cex.axis = 0.8, mgp = c(1,0.5,0), las = 1)
text(y=-1.5, x = b+0.3, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)
mtext(text="number of species", side=2, outer=F, line=2, las=0, cex=1)
mtext(text="Abundance class", side=1, outer=F, line=3.5)
mtext(text=c("b) Native focal species"), side=3, outer=F, line=5)

## 
#### Figure 2: trends in alpha richness effect size ####

par(mfrow = c(1,1), mar=c(0,0,2,1), oma=c(7,7,1,1))
sel = impsp.pos
M <- glmSRnat.overall
ylim=c(-50,200)

### Loop on selected species
for (i in 1:length(sel))  {
  # if (i ==4 | i == 7) plot.new()
  sp <- sel[i]
  es <- as.numeric(M$est[sp,])
  es <- exp(as.numeric(M$est[sp,]))*100-100
  n <- as.numeric(M$n.obs[sp,])[2:6]
  
  low <- as.numeric(M$CIlow[sp,])
  hi <- as.numeric(M$CIhi[sp,])
  low <- exp(as.numeric(M$CIlow[sp,]))*100-100
  hi <- exp(as.numeric(M$CIhi[sp,]))*100-100
  
  ## correct for small sample sizes
  small <- es
  small[which(n>=5)] <- NA
  es[which(n<5)] <- NA
  low[n<5] <- NA
  hi[n<5] <- NA
  
  # if ( !all(is.na(c(hi, low)))) lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
  # if ( all(is.na(c(hi, low)))) lims <-  max(c(abs(es)), na.rm=T) +0.02
  
  col <- rep("white",6)
  col[(1:6) >= as.numeric( M$impact.spread[sp,threshold])] <- "black"  # above threshold
  col[(2:6)[low<0]] <- NA  # robust negative coef
  col[(2:6)[n<5]] <- NA  # sufficient data points
  
  # plot background
  plot(1:5,  rep(0,5), ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F, bty = "n")
  
  #dotted line
  # abline(h=0,lty="dotted")
  abline(h=0,lty="dotted")
  
  # draw small sample sizes
  # if (!all(is.na(small))) points(1:5, small, pch=21, cex = 0.6)
  
  ## Add bootstrapped CI
  if ( !all(is.na(c(hi, low)))) {
    arrows(1:5,low,   1:5,hi, lwd=1, code =3, length=0.05, angle=90)
  }
  
  # draw the points and lines
  par(new=T)
  plot(c(1:5)[!is.na(es)],es[!is.na(es)], bg = col[2:6][!is.na(es)], pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F, bty = "n")
  
  #threshold line
  abline(v = M$impact.spread[sp,threshold]-1, lty="dashed")
  # th <-M$impact.spread[sp,threshold]
  # arrows(x0 =th-1, y0 = 2 ,x1 =th-1, y1 =1, length = 0.07, col = "grey60", lwd = 2, code = 2)
  
  # Add species name
  mtext(3, text=paste( species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.1, cex=0.7, line=0.2, las = 1)
  
  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  if ( i %in% c(1,2)) text(y=-75, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
  
  # Y axis
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, labels = F, lwd = 0, lwd.ticks = 1)
  if ( i %in% c(1,3)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  
  # Y axis name
  if ( i %in% c(1,3)) {
    mtext(2, text=substitute("% SRnat"[a], list(a = "[rare]")), adj=0.5, cex=0.8, line=2, las = 0)
  }
  
  box(bty = "o", lwd = 1)
  
}

mtext(1, text=c("Abundance class"), adj=0.5, line=5, las = 1, outer=T)
mtext(2, text=expression(paste("Effect size on native ", alpha,"-richness")), adj=0.5, line=4, las = 0, outer=T)


#### Figure 3: trends in gamma richness effect size ####

 par(mfrow = c(1,1), mar=c(0,0,2,1), oma=c(6,5,1,1))
ylim = c(-45,35)
xlim = c(0.5, 5.5)

for (i in 1 : length(sel)){
  
  # if (i ==4) plot.new()
  sp <- sel[i]
 
   n <- as.numeric(glmSRnat.overall$n.obs[sp,])[2:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)
  
  # create y axis with the gamma richness Standardized effect size
  y = (gamma.trend.nat$obs[sp,] - gamma.trend.nat$mean.null[sp,])
  small <- y
  small[which(n>=5)] <- NA
  y [ (2:6)[ n < 5] ]=NA
  y[1]=NA
  
  # plot background
  plot(x,y,ylim=ylim , xlim = c(0.5, 5.5), las = 1,type= "n", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols, bty="n")
  
  #  null expectation
  n <- as.numeric(glmSRnat.overall$n.obs[sp,])[1:6]
  lownull= as.numeric(gamma.trend.nat$q025.null[sp,]- gamma.trend.nat$mean.null[sp,] )
  lownull [ (1:6)[ n < 5] ]=NA  # remove classes with insufficient observations
  lownull <- lownull[2:6]       # remove the first abundance class for this graph
  hinull =  as.numeric(gamma.trend.nat$q975.null[sp,]- gamma.trend.nat$mean.null[sp,])
  hinull [ (1:6)[ n < 5] ]=NA   # remove classes with insufficient observations
  hinull <- hinull [2:6]        # remove the first abundance class for this graph
  
  #plot NULL CI permute.rare
  cinull <-c(lownull, hinull[5:1])
  cix <- c(x[2:6],x[6:2])
  cix[is.na(cinull)] <-NA
  polygon(na.omit(cix), na.omit(cinull), col="grey95", border="grey95")
  lines( c(1:5)[!is.na(lownull)], lownull[!is.na(lownull)],col="grey60", lty="dashed")
  lines( c(1:5)[!is.na(hinull)], hinull[!is.na(hinull)],col="grey60", lty="dashed")
  
  # if (!all(is.na(small))) points(1:6, small, pch=21, cex = 0.6)
  
  # annotations
  
  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  # if (i %in% c(8:11)) text(y=-53.2, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
  if (i %in% c(1,2)) text(y=-53.2, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
  
  # Y axis
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, labels = F, lwd = 0, lwd.ticks = 1)
  # if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  if ( i %in% c(1,3)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  
  # dotted horizontal
  abline(h=0,lty="dotted")
  
  # identify significant classes in black
  cols <- c(NA,  c(NA, "black") [  ((y[2:6] - hinull) >0) +1])
  
  
  #plot observed points and lines
  par(new=T)
  plot(x[!is.na(y)],y[!is.na(y)], bg = cols[!is.na(y)], pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F, bty="n")
  box(bty="o", lwd = 1)
  
  # Add arrow indicating threshold for alpha richness
  th <- glmSRnat.overall$impact.spread[sp,threshold]
  # arrows(x0 =th-1, y0 = as.numeric(y[th] +20) ,x1 =th-1, y1 = as.numeric(y[th] +10), length = 0.07, col = "grey60", lwd = 2, code = 2)
  arrows(x0 =th-1, y0 = -35,x1 =th-1, y1 = -42, length = 0.07, col = "grey60", lwd = 1, code = 1)
  
  # # Add spearman test if more than 2 points :
  # if ( length(na.omit(as.numeric(y)))>2) {
  # fit <- cor.test(2:6,as.numeric(y[2:6]), method = "spearman", exact = FALSE) 
  # mtext(3, text=paste(round(fit$estimate,2), p2star(fit$p.value)) ,
  #       font = 1, outer= F,adj=0.99, cex=0.7, line=-1, las = 1)
  # }
  
  # Add species name
  mtext(3, text=paste( species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.1, cex=0.7, line=0.2, las = 1)
}

mtext(1, text=c("Abundance class"), adj=0.5, line=5, las = 1, outer=T)
mtext(2, text=expression(paste("Effect size on native ", gamma,"-richness")), adj=0.5, line=3, las = 0, outer=T)




