# Figure Article thresholds of impact

abclasses= c("Occasional" ,"Common-Occasional",  "Common", "Abundant-Common", "Abundant","Dominant")
abclasses= c("Rare" ,"Occasional",  "Frequent", "Common", "Abundant","Dominant")

############## Barplot of frequencies #############

effects =list(glmSRnat.sum$class.summary,glmSRali.sum$class.summary )     

n <- length(unique(effects[[1]]$group)) 

par(mfcol=c(n,2), oma=c(5,6,4,1), mar=c(1,1,1,2), las=1)
for (j in 1:2) {
  sum.df = effects[[j]] 
  for (i in 1:n) {
    S <- as.data.frame(sum.df[sum.df$group == unique(sum.df$group)[i],])
    S$prop.negative.target <- S$n.negative.target/S$n.target
    
    barplot(S$nb.sp, ylim=c(0,max(40, S$nb.sp)),
            col= "grey",  border= NA, axes=F)
    par(new=T)
    b <- barplot(S$freq.thr, col="black" , ylim=c(0,max(40, S$nb.sp)), border= NA, axes=F)
    axis(2, tcl = -0.3, cex.axis = 0.8, mgp = c(1,0.5,0))
    if (i==2) text(y=-3, x = b, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)
    
    if (j==2 & i==1) {
      legend(x=0.3, y=40, bty="n", bg="white",legend=c("Critical abundance", "Total occurrences"),
             fill=c("black",  "grey"), border= c("black",  "grey"), cex=0.9, xpd=NA)
    }
    if(j==1) {
      mtext(text="number of species", side=2, outer=F, line=1.5, las=0, cex=0.85) 
    }
  }   
}

mtext(side=2, text=c("Alien\nfocal\nspecies", "Native\nfocal\nspecies"),line=4, at=c(0.25,0.75),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("a) Native richness", "b) Alien richness"),line=0, at=c(0.25,0.75),adj=0.5, outer=T, las=1)

mtext(text="Abundance class", side=1, outer=T, line=3)

    
#########  Details per significant species ##############
threshold = "th.CI"

# Effect size for significant species

# x11()
# plot.glm(M = glmSRnat.overall, var= "SRnat",  db= db,sel.criteria = "sp.target", type= "overall.boot",
#          ES = T, threshold = threshold, sp = sp.target,  panel = c(3,4), ylim=c(-2,2))
# mtext(2,text = "effect size on Native richness", outer=T, line= 1)
# mtext(1,text = "abundance class", outer=T, line= 1)
# 
# x11()
plot.glm(M = glmSRali.overall, var= "SRali",  db= db,sel.criteria = "sp.target", type= "overall.boot",
         ES = T, threshold = threshold, sp = sp.target, panel = c(3,4), ylim=c(-0.6,0.6))
mtext(2,text = "effect size on Alien richness", outer=T, line= 1)
mtext(1,text = "abundance class", outer=T, line= 1)



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
    mtext(2, text="GLM coefficient", ,adj=0.5, cex=0.7, line=1.5, las = 0)
    }
   }

mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("Effect on native alpha richness"), adj=0.5, line=4, las = 0, outer=T)

##### gamma trends   ############

par(mfrow = c(3,4), mar=c(0,0,1,1), oma=c(6,6,2,1))
ylim=c(-5,2)
for (i in 1 : length(impsp)){
  
  if (i ==4) plot.new()
  sp <- impsp[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ (gamma.trend.nat$P.gamma[i,] <=0.025) +1]
  n <- as.numeric(M$n.obs[sp,])[2:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)
  
  # create y axis with the gamma richness Standardized effect size
  y = mnnd.above$z[1,]
  y [ (2:6)[ n < 5] ]=NA
  y[1]=NA
  
  
  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8,9,10,11)) text(y=-5.2, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)

  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)  
  abline(h=0,lty="dotted")
  
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


##### Gamma and alpha on same graph ########

par(mfrow = c(3,4), mar=c(0,0,2,1), oma=c(6,6,2,3))
M <- glmSRnat.overall
ylim=c(-3,2)
### Loop on selected species 
for (i in 1:length(sel))  {
  if (i == 4) {
    plot(1:10, 1:10, ann=F, axes = F, type ="n")
  legend(3,8, legend =c(expression(alpha*"-richness"), expression(gamma*"-richness")),
         pch =c(21, 24) , col =c("black", "grey50"),pt.bg =c("black", "grey50"), lty = "solid", bty = "n")
  }
  sp <- sel[i]
  es <- as.numeric(M$est[sp,])
  n <- as.numeric(M$n.obs[sp,])[2:6]
  
  low <- as.numeric(M$CIlow[sp,])
  hi <- as.numeric(M$CIhi[sp,])
  M$thresh <-M$impact.spread
  
  ## correct for small sample sizes
  small <- es
  small[which(n>=5)] <- NA
  es[which(n<5)] <- NA
  low[n<5] <- NA
  hi[n<5] <- NA
  
  if ( !all(is.na(c(hi, low)))) lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
  if ( all(is.na(c(hi, low)))) lims <-  max(c(abs(es)), na.rm=T) +0.02
  
  col <- rep("white",6)
  col[(1:6)>= as.numeric( M$thresh[sp,threshold])] <- "black"  # above threshold
  col[(2:6)[hi>0]] <- "white" # robust negative coef
  col[(2:6)[n<5]] <- "white" # sufficient data points
  
  ### prepare background plots
  plot(1:5,  rep(0,5),ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F)
  
  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(8:11)) text(y=-3.4, x = 1:5, labels= abclasses[2:6],  cex=1, srt=50, adj=1, xpd = NA)
  axis(1, tcl= 0.2,  mgp=c(1,0.5,0), las=1, labels = FALSE) 
  
  if ( i %in% c(1,4,8)) axis(2, tcl= 0,   mgp=c(1,0.5,0), las=1, cex.axis = 0.9, adj =1)  
   axis(2, tcl= 0.2, label = FALSE ,  mgp=c(1,0.5,0)) 
  
  abline(h=0,lty="dotted")
  abline(v = M$thresh[sp,threshold]-1, lty="dashed")
  
  if (!all(is.na(small))) points(1:5, small, pch=".")  
  
  ## Add bootstrapped CI
  if ( !all(is.na(c(hi, low)))) {
    arrows(1:5,low,   1:5,hi, lwd=1, code =3, length=0.02, angle=90)
  }
  
  # draw the points and lines
  par(new=T)
  plot(1:5,es, bg = col[2:6], pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), 
       lwd = 1.5, cex = 1.2,
       bty="n", xaxt = "n", yaxt="n", ann=F)
  
  
  ### gamma
  # identify significant classes in black
  bgs <- c(NA, "grey60") [ (gamma.trend.nat$P.gamma[i,] <=0.025) +1]
  cols <- c("grey60", "grey60") [ (gamma.trend.nat$P.gamma[i,] <=0.025) +1]
  
  # create x axis for the 5 abundance classes
  x = c(0:5)
  
  # create y axis with the gamma richness Standardized effect size
# #   y =(gamma.trend.nat$obs[i,] - gamma.trend.nat$mean[i,])/ gamma.trend.nat$sd[i,]
#   y =  (gamma.trend.nat$obs[i,] - gamma.trend.nat$mean[i,])/ gamma.trend.nat$mean[i,]
y =  (gamma.trend.nat$obs[i,] - gamma.trend.nat$mean[i,])
  y [ (2:6)[ n < 5] ]=NA
  y[1]=NA
  par(new=T)
#   plot(x,y, bg = cols,col =cols, pch=24, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)
plot(x,y, bg = bgs,col=cols, pch=24, type = "b", cex = 1.2,
     ylim = c(-45,30), xlim = c(0.5, 5.5), 
     bty="n", xaxt = "n", yaxt="n", ann=F)
if ( i %in% c(3,7,11)) axis(4, tcl= 0,  mgp=c(1,0.3,0), las=1, cex.axis = 0.9, adj =1) 
axis(4, tcl= 0.2,  mgp=c(1,0.5,0), labels = FALSE) 

  
  # Add species name
  mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)
  
  # Y axis label
  if ( i %in% c(1,4,8)) {  
    mtext(2, text=expression("ES("*alpha*")"), ,adj=0.5, cex=0.9, line=2, las = 0)
  }
  if ( i %in% c(3, 7,11)) {  
    mtext(4, text=expression("ES("*gamma*")"), ,adj=0.5, cex=0.9, line=2, las = 0)
  }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4.5, las = 1, outer=T)
mtext(2, text=c("Effect size on native richness"), adj=0.5, line=4.1, las = 0, outer=T)

######### ILLUSTRATION of 3 species alpha and gamma richness   ############
threshold ="th.CI"
sp.target = c("ACHMIL","LOLPER", "ANTODO")
sel <-  sp.target

x11()

par(mfrow = c(2,3), mar=c(2,2,1,1), oma=c(2,9,2,1))

M=glmSRnat.overall
ylim=c(-3,3)
### Loop on selected species 
for (sp in sel)  {
  
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
  col[ as.numeric( M$thresh[sp,threshold])] <- "black"  # identify the threshold
  
  
  plot(1:6,  rep(0,6),ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F)
  axis(1,at = 1:5, label = names(M$est), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)  
  abline(h=0,lty="dotted")
  
  # if (!all(is.na(small))) points(1:5, small, pch=".")  
  
  ## Add bootstrapped CI
  if ( !all(is.na(c(hi, low)))) {
    arrows(1:5,low,   1:5,hi, lwd=1, code =3, length=0.05, angle=90,  xlim = c(0.5, 5.5),ylim=ylim)
  }
  par(new=T)
  plot(1:5,es, bg = col[2:6], pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)
  
  # mtext(side=3, line=-1.2, at=1:5,text=n, cex=0.7, las=1, col = "darkgrey")
  
  mtext(3, text=paste(species[sp, "Genus"], " ",species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.5, cex=0.8, line=0.5, las = 1)
  if(sp == "ACHMIL")   mtext(2, text="GLM coefficient", adj=0.5, cex=0.7, line=2, las = 0)
}

for (i in sel){
  cols <- c("white", "black") [ (gamma.trend.nat$P.gamma[i,] <=0.025) +1]
  x = 1:6
  y =(gamma.trend.nat$obs[i,] - gamma.trend.nat$mean[i,])/ gamma.trend.nat$sd[i,]
  y[1]=NA
  plot(x,y,ylim= c(-5,4),xlim= c(1,6.3),ann=F, las = 1, pch = 21, col ="black",bg = cols, type= "b")
  abline(h=0,lty="dotted")
  if(i == "ACHMIL")   mtext(2, text="SES",adj=0.5, cex=0.7, line=2, las = 0)
  
}

mtext(2, at = c(0.77, 0.3),text=c("Alpha\nrichness", "Gamma\n richness"),
      adj=0.5,  line=5, las = 1, outer=T)
mtext(1, text=c("Abundance class"), adj=0.5, line=1, las = 1, outer=T)

# ### frequency with bootstrap variance #######
x11()

tab <- apply(glmSR.overall$crit.vals[,], 2, FUN=function(x) {
  x <- factor(x,levels= c("2","3","4","5","6"))
  f = table(x)
})

sd = apply(tab,1,  FUN= sd, na.rm=T)
barplot(tab[match(2:6, rownames(tab)),1], ylim= c(0,max(tab) +2))
segments(1:4, tab[,1] - sd, 1:4, tab[,1] + sd)
mtext(2,text = "Frequency critical value", outer=F, line= 3)
mtext(1,text = "abundance class", outer=F, line= 3)


##### gamma above vs. below trends   ############

par(mfrow = c(3,4), mar=c(0,0,1,1), oma=c(6,6,2,1))

list.data <-divpart.nat.perm
gamma.loss <- as.data.frame(matrix(NA, nrow = length(impsp), ncol =6))
rownames(gamma.loss) <- impsp

ylim=c(-1,1)
for (i in 1 : length(impsp)){
  
  if (i ==4) plot.new()
  sp <- impsp[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ (list.data$P.above[i,] <=0.025) +1]
  n <- as.numeric(M$n.obs[sp,])[2:6]
  # create x axis for the 5 abundance classes
  x = c(0:5)
  
  # create y axis with the gamma richness Standardized effect size
  y =  (list.data$gamma.above[i,] - list.data$null.above[i,])/ list.data$null.above[i,]
  y [ (2:6)[ n < 5] ]=NA
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
mtext(2, text=c("Effect size of native gamma richness"), adj=0.5, line=3, las = 0, outer=T)

######### Graphs of delta gamma vs delta alpha  ###########
#community change plots
## species order
cm <- comm[(rownames(comm) %in% realgrasslands),]
cm <- cm[,colnames(cm)%in%natives]
cm <- ceiling(cm>0)
cm <- cm[,colSums(cm)>0]
ord <- colnames(cm)[order(colSums(cm))]
ord <- ord[ord %in% natives]

## NON STANDARDIZED effect sizes
dg<- divpart.nat.perm$obs$above.gamm - divpart.nat.perm$null$above.gamm
da<- divpart.nat.perm$obs$above.alpha - divpart.nat.perm$null$above.alpha
dga<- divpart.ali.perm$obs$above.gamm - divpart.ali.perm$null$above.gamm
daa<- divpart.ali.perm$obs$above.alpha - divpart.ali.perm$null$above.alpha

cgam <- c("grey50", "black")[ (divpart.nat.perm$P$above.gamm<=0.05 | divpart.nat.perm$P$above.gamm>=0.95) + 1]

### figure 4 panels : ##########
x11()

dg<- -(table.div.part$GRc - table.div.part$GRnull )
# dg<- -(table.div.part$GRc - table.div.part$GRnull )/table.div.part$GRsdnull
da<-  -( table.div.part $aRc -table.div.part$aRo )
cgam <- c("grey50", "black")[ (table.div.part$GR.P<=0.05) + 1]

par(mfrow = c(2,2), mar=c(4,4,2,2), las = 1)
plot(da, dg ,pch = 20, col=cgam , ann=F,ylim=c(1,40), xlim= c(1,10), type ="n")

text(x = da,y = dg,label = rownames(divpart.nat.perm$P) , col=cgam, cex =0.85)

f <- cor.test(da ,dg)
mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.8)

mtext(1,text = expression(Delta*" alpha richness"), line=2.5, cex =0.85 )
mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)
mtext(3,text = "a)", line=0.5  , las=0, adj = 0, cex =0.85)


## threshold
thresh <-  glmSRnat.overall$impact.spread[impsp,  "th.CI"]
plot (thresh ,dg , pch = 20, col=cgam , ann=F, xlim= c(1,6.5), type = "n")
f <- cor.test(thresh ,dg)

text(x = thresh , y = dg,
     label = rownames(divpart.nat.perm$P) ,col=cgam, cex =0.85 )

mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.8)

mtext(1,text = "Critical abundance", line=2.5, cex =0.85 )
mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)
mtext(3,text = "b)", line=0.5  , las=0, adj = 0, cex =0.85)

### gamma above/below against threhold values
spread <-  glmSRnat.overall$impact.spread[impsp,  "n.plot.impact"]
plot (spread,dg , pch = 20, col=cgam , ann=F, log = "x", xlim =c(2.5, 700), type= "n")
f <- cor.test(spread ,dg)
text(x = spread,
     y = dg,
     label = rownames(divpart.nat.perm$P) ,col=cgam, cex =0.85)

mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.8)

mtext(1,text = "Number of plots > critical abundance", line=2.5, cex =0.85 )
mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)
mtext(3,text = "c)", line=0.5  , las=0, adj = 0, cex =0.85)

# spatial clustering
cluster<-  mnnd.th[impsp, "P<obs_above"]
cluster<- - (mnnd.th[impsp, "mnnd.obs_above"] - mnnd.th[impsp, "null.mean_above"])/ mnnd.th[impsp, "null.sd_above"] 
plot (cluster,dg , pch = 20, col=cgam , ann=F,  type="n",xlim=c(-1,2.6))
f <- cor.test(cluster ,dg)

text(x = cluster,
     y = dg,
     label = rownames(divpart.nat.perm$P) ,col=cgam, cex =0.7 )

mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.8)

mtext(1,text = "Spatial clustering of plots > critical abundance", line=2.5, cex =0.85 )
mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)
mtext(3,text = "d)", line=0.5  , las=0, adj = 0, cex =0.85)

######### alternative figures : ###########




x11()
# prop of impact
spread <-  glmSRnat.overall$impact.spread[impsp,  "prop.plot.impact"]
plot (spread,dg , pch = 20, col=cgam , ann=F, xlim =c(0, 1), type="n")
f <- cor.test(spread ,dg)

text(x = spread,
     y = dg,
     label = rownames(divpart.nat.perm$P) ,col=cgam, cex =0.7 )

mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 0, cex=0.8)

mtext(1,text = "Prop. of plots > critical abundance", line=2.5, cex =0.85 )
mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0, cex =0.85)


x11()
###  clustering vs threshold
cluster<-  (mnnd.th[impsp, "mnnd.obs_above"] - mnnd.th[impsp, "null.mean_above"])/ mnnd.th[impsp, "null.sd_above"] 
thresh <-  glmSRnat.overall$impact.spread[impsp,  "th.CI"]

plot (thresh,cluster, pch = 20, col=cgam , ann=F,  type="n")
f <- cor.test(cluster ,thresh, method = "spearman")

text(x = thresh,
     y = cluster,
     label = rownames(divpart.nat.perm$P) ,col=cgam, cex =0.7 )

mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 0, cex=0.8)

mtext(1,text = "critical abundance", line=2.5, cex =0.85 )
mtext(2,text ="SES of MNND", line=2.5  , las=0, cex =0.85)

x11()
###  alpha loss vs threshold
da<-  -( table.div.part $aRc - table.div.part$aRo )
da<-  as.numeric(impact.SRnat[impsp,]$wtd.mean.dif)
da<-  as.numeric(impact.SRnat[impsp,]$th.dif)

thresh <-  glmSRnat.overall$impact.spread[impsp,  "th.CI"]

plot (thresh,da, pch = 20, col=cgam , ann=F,  type="n")
f <- cor.test(da,thresh, method = "spearman")

text(x = thresh,
     y = da,
     label = impsp ,col="black", cex =0.7 )

mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 0, cex=0.8)

mtext(1,text = "critical abundance", line=2.5, cex =0.85 )
mtext(2,text ="alpha loss", line=2.5  , las=0, cex =0.85)


### spatial clustering vs loss in gamma richness at each abundance  ############


par(mfrow = c(3,4), mar=c(0,0,1,1), oma=c(6,6,2,1))

list.data <-divpart.nat.perm
gamma.loss <- as.data.frame(matrix(NA, nrow = length(impsp), ncol =6))
rownames(gamma.loss) <- impsp

ylim=c(-50,20)
for (i in 1 : length(impsp)){
  
  if (i ==4) plot.new()
  sp <- impsp[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ (list.data$P.above[i,] <=0.025) +1]
  cols [2] <- "red"
  n <- as.numeric(M$n.obs[sp,])[2:6]
  # create x axis for the 5 abundance classes

  # create y axis with the gamma richness Standardized effect size
  
  x= mnnd.above$z[i,]
  
  y =  (list.data$gamma.above[i,] - list.data$null.above[i,])
  y [ (2:6)[ n < 5] ]=NA
  y[1]=NA
  
  gamma.loss[sp,] = y
  
  # plot background
  plot(x,y,ylim=ylim, xlim = c(-5,2), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)
  
  axis(1, tcl= 0.1,mgp=c(1,0.5,0),las=1)
    
  if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)  
  abline(h=0,lty="dotted")
  abline(v=0,lty="dotted")
  
  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(-5,2), xaxt = "n", yaxt="n", ann=F)
  
  # Add species name
  mtext(3, text=paste(species[sp, "Genus"], "\n",species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=-2.5, las = 1)
  
  #   # Y axis label
  #   if ( i %in% c(1,4,7)) {  
  #     mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
  #   }
}

mtext(1, text=c("SES MNND"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("Effect size of native gamma richness"), adj=0.5, line=3, las = 0, outer=T)

