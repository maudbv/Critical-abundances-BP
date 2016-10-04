# Define some elements
abclasses= c("Rare" ,"Occasional",  "Frequent", "Common", "Abundant","Dominant")
threshold ="th.CI"
sel <- impsp
sel <- sel[c(1,10,11,2:9)]


### Figure 1: barplot of frequencies:    #########

# extract frequency table:
effects =list(SRnat = glmSRnat.sum$class.summary,SRali = glmSRali.sum$class.summary )
n <- length(unique(effects[[1]]$group))
for (j in 1:2) {
  sum.df = effects[[j]]
  for (i in 1:n) {
    S <- as.data.frame(sum.df[sum.df$group == unique(sum.df$group)[i],])
    }
}

#draw barplot
par(mar=c(5,3,2,2), las=1)
S <- effects$SRnat[effects$SRnat$group == "ALIEN:1",]
barplot(S$nb.sp, ylim=c(0,max(40, S$nb.sp)),col= "grey80",  border= NA, axes=F)
par(new=T)
barplot(S$freq.negative.above, ylim=c(0,max(40, S$nb.sp)),col= "grey60",  border= NA, axes=F)
par(new=T)
b <- barplot(S$freq.thr, col="black" , ylim=c(0,max(40, S$nb.sp)), border= NA, axes=F)
axis(2, tcl = -0.3, cex.axis = 0.8, mgp = c(1,0.5,0), las = 1)
text(y=-1.5, x = b, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)
legend(x=3.6, y=48, bty="n", bg="white",legend=c( "All species", "with negative effects","at critical abundance"),
             fill=c( "grey80","grey60","black"), border= c("grey90","grey60","black"), cex=0.85, xpd=NA,y.intersp =1.5)
mtext(text="number of target species", side=2, outer=F, line=2, las=0, cex=1)
mtext(text="Abundance class", side=1, outer=F, line=3.5)

#### Figure 2: trends in alpha richness effect size ####

par(mfrow = c(3,4), mar=c(0,0,2,2), oma=c(7,6,2,3))

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
  
  ## correct for small sample sizes
  small <- es
  small[which(n>=3)] <- NA
  es[which(n<3)] <- NA
  low[n<3] <- NA
  hi[n<3] <- NA
  
  if ( !all(is.na(c(hi, low)))) lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
  if ( all(is.na(c(hi, low)))) lims <-  max(c(abs(es)), na.rm=T) +0.02
  
  col <- rep("white",6)
  col[(1:6)>= as.numeric( M$impact.spread[sp,threshold])] <- "black"  # above threshold
  col[(2:6)[hi>0]] <- NA  # robust negative coef
  col[(2:6)[n<5]] <- NA  # sufficient data points
  
  # plot background
  plot(1:5,  rep(0,5), ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F, bty = "n")
  
  #dotted line
  abline(h=0,lty="dotted")
  
 
  # draw small sample sizes
  if (!all(is.na(small))) points(1:5, small, pch=".")
  
  ## Add bootstrapped CI
  if ( !all(is.na(c(hi, low)))) {
    arrows(1:5,low,   1:5,hi, lwd=1, code =3, length=0.05, angle=90)
  }
  
  # draw the points and lines
  par(new=T)
  plot(1:5,es, bg = col[2:6], pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F, bty = "n")
  
  #threshold line
  abline(v = M$impact.spread[sp,threshold]-1, lty="dashed")
  # th <-M$impact.spread[sp,threshold]
  # arrows(x0 =th-1, y0 = 2 ,x1 =th-1, y1 =1, length = 0.07, col = "grey60", lwd = 2, code = 2)
  
  # Add species name
  mtext(3, text=paste(letters[i],") ", species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.1, cex=0.8, line=0.2, las = 1)
  
  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  if (i %in% c(8:11)) text(y=-3.6, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
  
  # Y axis
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, labels = F, lwd = 0, lwd.ticks = 1)
  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)

   # Y axis name
  if ( i %in% c(1,4,8)) {
    mtext(2, text="GLM coefficient", adj=0.5, cex=0.8, line=1.5, las = 0)
  }
  
  box(bty = "o", lwd = 1)
  
}

mtext(1, text=c("Abundance class"), adj=0.5, line=5, las = 1, outer=T)
mtext(2, text=expression(paste("Effect size for native ", alpha,"-richness")), adj=0.5, line=4, las = 0, outer=T)

#### Figure 3: trends in gamma richness effect size ####

par(mfrow = c(3,4), mar=c(0,0,2,2), oma=c(7,6,2,3))
ylim = c(-45,35)
xlim = c(0.5, 5.5)

for (i in 1 : length(sel)){
  
  if (i ==4) plot.new()
  sp <- sel[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ (gamma.trend.nat$P.gamma[sp,] <=0.025) +1]
  n <- as.numeric(glmSRnat.overall$n.obs[sp,])[2:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)
  
  # create y axis with the gamma richness Standardized effect size
  y = (gamma.trend.nat$obs[sp,] - gamma.trend.nat$mean.null[sp,])
  y [ (2:6)[ n < 5] ]=NA
  y[1]=NA
  
  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "n", xaxt = "n", yaxt="n", ann=F,
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
  lines(1:5, lownull,col="grey60", lty="dashed")
  lines(1:5, hinull,col="grey60", lty="dashed")
  
  # annotations

  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  if (i %in% c(8:11)) text(y=-53.2, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
  
  # Y axis
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, labels = F, lwd = 0, lwd.ticks = 1)
  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  
  # dotted horizontal
  abline(h=0,lty="dotted")
  
  #plot observed points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F, bty="n")
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
  mtext(3, text=paste(letters[i],") ", species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.1, cex=0.8, line=0.2, las = 1)
}

mtext(1, text=c("Abundance class"), adj=0.5, line=5, las = 1, outer=T)
mtext(2, text=expression(paste("Effect size for native ", gamma,"-richness")), adj=0.5, line=3, las = 0, outer=T)


### Figure 4: correlation of deltagamma with delta alpha and with spread  ######

# loss in gamma native richness
dg<- table.div.part$deltagamma - table.div.part$delta.null.permute.all
dg1<- table.div.part$deltagamma - table.div.part$delta.null.beta
dg2<- table.div.part$GRc - table.div.part$GRnull

# loss in alpha native richness
da<- table.div.part$aRc -  table.div.part$aRo
da1<- table.div.part$aRc -  table.div.part$aRnull

#spatial distributions:
spread <- table.div.part$n.plot.impact
prevalence <- table.div.part$prevalence

##colour code
# cgam <- c("grey50", "black")[ (table.div.part$GRP.permute.all<=0.05) + 1]
cgam <- "black"

#loss in alpha vs. loss in gamma
par(mfrow = c(1,2), mar=c(3,1,2,1), oma=c(1,2,0,0), las = 1)
plot(da, dg ,pch = 20, col=cgam , ann=F, xlim= c(-4,0), type ="n", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x = da,y = dg,label = rownames(table.div.part) , col=cgam, cex =0.7)
f <- cor.test(da ,dg1)
box(bty="o")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8, font = 3)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1)
mtext(1,text = expression(Delta*alpha*"-richness"), line=2 )
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=2 , las=0)

### spread above critical abundance
print(f <- cor.test(spread ,dg))
plot (spread,dg , pch = 20, col=cgam , ann=F, log = "x", xlim =c(2.5, 700), type ="n", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=F, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x = spread, y = dg,label = rownames(table.div.part),col=cgam, cex=0.7)
box(bty="o")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8)
mtext(1,text = "Number of plots > critical abundance", line=2 )
mtext(3, text = 'b)',adj = 0, cex=0.8, font = 1)

##Extra figure #####

### gamma vs. prevalence
par(mfrow = c(1,1), mar=c(4,3,2,1), oma=c(0,0,0,0), las = 1)
print(f <- cor.test(prevalence,dg))
plot (prevalence,dg , pch = 20, ann=F, log = "x", xlim =c(100, 900), type= "n")
text(x = prevalence, y = dg,label = rownames(table.div.part),  col=cgam ,cex=0.7)
mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.8)
mtext(1,text = "Number of plots", line=2)
mtext(2,text = expression(Delta*" gamma richness"), line=2  , las=0)
