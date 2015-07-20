# Figure Article thresholds of impact

abclasses= c("Occasional" ,"Common-Occasional",  "Common", "Abundant-Common", "Abundant","Dominant")
abclasses= c("Rare" ,"Occasional",  "Frequent", "Common", "Abundant","Dominant")

############## Barplot of frequencies #############

effects =list(glmSRnat.sum$class.summary,glmSRali.sum$class.summary )     

n <- length(unique(effects[[1]]$group)) 

par(mfcol=c(n,2), oma=c(5,6,4,1), mar=c(1,2,1,2), las=1)
for (j in 1:2) {
  sum.df = effects[[j]] 
  for (i in 1:n) {
    S <- as.data.frame(sum.df[sum.df$group == unique(sum.df$group)[i],])
    S$prop.negative.target <- S$n.negative.target/S$n.target
    
    barplot(S$nb.sp, ylim=c(0,max(40, S$nb.sp)),
            col= "grey",  border= NA)
    par(new=T)
    b <- barplot(S$freq.thr, col="black" , ylim=c(0,max(40, S$nb.sp)), border= NA)
    
    if (i==2) text(y=-3, x = b, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)
    
    if (j==2 & i==1) {
      legend(x=0.3, y=40, bty="n", bg="white",legend=c("Critical abundance", "total observed"),
             fill=c("black",  "grey"), border= c("black",  "grey"), cex=0.8, xpd=NA)
    }
    if(j==1) {
      mtext(text="number of species", side=2, outer=F, line=2.5, las=0, cex=0.85) 
    }
  }   
}

mtext(side=2, text=c("Alien\nfocal\nspecies", "Native\nfocal\nspecies"),line=4, at=c(0.25,0.75),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("Native richness", "Alien richness"),line=0, at=c(0.25,0.75),adj=0.5, outer=T, las=1)
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

x11()
par(mfrow = c(3,3), mar=c(0,0,1,1), oma=c(6,6,2,1))
M <- glmSRnat.overall
  ylim=c(-3,3)
  ### Loop on selected species 
  for (i in 1:length(sel))  {
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
    col[ as.numeric( M$thresh[sp,threshold])] <- "black"  # identify the threshold
    
    
    plot(1:5,  rep(0,5),ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F)

    axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
    if (i %in% c(7,8,9)) text(y=-3.6, x = 1:5, labels= abclasses[2:6],  cex=0.9, srt=45, adj=1, xpd = NA)
    
    if ( i %in% c(1,4,7)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)  
    abline(h=0,lty="dotted")
    
    if (!all(is.na(small))) points(1:5, small, pch=".")  
    
    ## Add bootstrapped CI
    if ( !all(is.na(c(hi, low)))) {
      arrows(1:5,low,   1:5,hi, lwd=1, code =3, length=0.05, angle=90)
    }
    
    # draw the points and lines
    par(new=T)
    plot(1:5,es, bg = col[2:6], pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)
    
    # Add species name
    mtext(3, text=paste(species[sp, "Genus"], "\n",species[sp, "Species"], sep="") ,
                       font = 3, outer= F,adj=0.9, cex=0.7, line=-2.5, las = 1)
    
    # Y axis label
    if ( i %in% c(1,4,7)) {  
    mtext(2, text="GLM coefficient", ,adj=0.5, cex=0.7, line=1.5, las = 0)
    }
   }

mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("Effect on native alpha richness"), adj=0.5, line=4, las = 0, outer=T)

### gamma trends

par(mfrow = c(3,3), mar=c(0,0,1,1), oma=c(6,6,2,1))
ylim=c(-5,4)
for (i in 1 : length(impsp)){
  sp <- impsp[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ (gamma.trend.nat$P.gamma[i,] <=0.025) +1]
  
  # create x axis for the 5 abundance classes
  x = c(0:5)
  
  # create y axis with the gamma richness Standardized effect size
  y =(gamma.trend.nat$obs[i,] - gamma.trend.nat$mean[i,])/ gamma.trend.nat$sd[i,]
  y[1]=NA
  
  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)

  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  if (i %in% c(7:9)) text(y=-6, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)

  if ( i %in% c(1,4,7)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)  
  abline(h=0,lty="dotted")
  
  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)
    
  # Add species name
  mtext(3, text=paste(species[sp, "Genus"], "\n",species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=-2.5, las = 1)
  
  # Y axis label
  if ( i %in% c(1,4,7)) {  
    mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
  }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
mtext(2, text=c("Variation in native gamma richness"), adj=0.5, line=4, las = 0, outer=T)



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



######### Graphs of delta gamma vs delta alpha  ###########
#community change plots
## species order
cm <- comm[(rownames(comm) %in% realgrasslands),]
cm <- cm[,colnames(cm)%in%natives]
cm <- ceiling(cm>0)
cm <- cm[,colSums(cm)>0]
ord <- colnames(cm)[order(colSums(cm))]
ord <- ord[ord %in% natives]



# gams <-  rowSums(divpart.nat.perm$obs[,c("shared","lost", "gained")])
# alphas <- as.numeric(impact.SRnat[rownames(divpart.nat.perm$obs),]$SRo)
# dg<- (divpart.nat.perm$obs$above.gamm - divpart.nat.perm$obs$below.gamm)/gams
# da<- (divpart.nat.perm$obs$above.alpha - divpart.nat.perm$obs$below.alpha) / alphas
# gams <-  rowSums(divpart.ali.perm$obs[,c("shared","lost", "gained")])
# alphas <- as.numeric(impact.SRali[rownames(divpart.nat.perm$obs),]$SRo)
# dga<- (divpart.ali.perm$obs$above.gamm - divpart.ali.perm$obs$below.gamm)/gams
# daa<- (divpart.ali.perm$obs$above.alpha - divpart.ali.perm$obs$below.alpha) / alphas

## NON STANDARDIZED effect sizes
dg<- divpart.nat.perm$obs$above.gamm - divpart.nat.perm$null$above.gamm
da<- divpart.nat.perm$obs$above.alpha - divpart.nat.perm$null$above.alpha
dga<- divpart.ali.perm$obs$above.gamm - divpart.ali.perm$null$above.gamm
daa<- divpart.ali.perm$obs$above.alpha - divpart.ali.perm$null$above.alpha


cgam <- c("grey50", "black")[ (divpart.nat.perm$P$above.gamm<=0.05 | divpart.nat.perm$P$above.gamm>=0.95) + 1]


x11()
par( mar=c(4,4,2,2), las = 1)
plot(da, dg ,xlim =c(-5,2),ylim =c(-40,7) , pch = 3, col=cgam , ann=F)

# segments (da, dg, daa, dga,col="grey", lty ="dotted")
# points(daa, dga, pch = 21, col="grey" ,bg = "grey"  )

text(x = da,
     y = dg-01,
     label = rownames(divpart.nat.perm$P) , col=cgam, cex =0.7 )

abline(h=0, v=0, col="grey")
# legend("topright", legend = c("native richness", "alien richness"), cex=0.85,
#        pch= c(21,3), pt.bg=c("black","white"), bty="o")

mtext(1,text = expression(Delta*" alpha richness"), line=2.5 )
mtext(2,text = expression(Delta*" gamma richness"), line=2.5  , las=0)

# PERCENTS effect sizes
# dg<- (divpart.nat.perm$obs$above.gamm - divpart.nat.perm$null$above.gamm)/divpart.nat.perm$null$above.gamm
# da<- (divpart.nat.perm$obs$above.alpha - divpart.nat.perm$null$above.alpha)/divpart.nat.perm$null$above.alpha
# dga<- (divpart.ali.perm$obs$above.gamm - divpart.ali.perm$null$above.gamm)/divpart.ali.perm$null$above.gamm
# daa<- (divpart.ali.perm$obs$above.alpha - divpart.ali.perm$null$above.alpha)/divpart.ali.perm$null$above.alpha
# 
# 
# x11()
# par( mar=c(4,4,2,2))
# plot(da, dg ,xlim =c(-1,1),ylim =c(-1,1) , pch = 21, col=cgam ,bg = cgam, ann=F)
# # points(daa, dga,  pch = 21, col=cgam )
# # segments (da, dg, daa, dga,col=cgam)
# 
# text(x = da,
#      y = dg-0.01,
#      label = rownames(divpart.nat.perm$P) , col=cgam, cex =0.7 )
# 
# abline(h=0, v=0, col="grey")
# legend("topright", legend = c("Native richness", "Alien Richness"), cex=0.9,
#        pch= c(21,21), pt.bg=c("black","white"), bty="n")
# 
# mtext(1,text = "Change in Alpha", line=2.5 )
# mtext(2,text = "Change in Gamma", line=2.5  )
# 

# 
## Effect sizes :   STANDARDIZED
dg<- divpart.nat.perm$z$above.gamm 
da<- divpart.nat.perm$z$above.alpha

dga<- divpart.ali.perm$z$above.gamm
daa<- divpart.ali.perm$z$above.alpha

cgam <- c("darkgrey", "black")[ (divpart.nat.perm$P$above.gamm<=0.05 | divpart.nat.perm$P$above.gamm>=0.95) + 1]
# 
x11()
par( mar=c(4,4,2,2), las = 1)
plot(daa, dga ,xlim =c(-8,2),ylim =c(-5,2) , pch = 3, col="darkgrey", ann=F)
segments (da, dg, daa, dga,col="grey", lty ="dotted")
points(da, dg, pch = 21, col=cgam ,bg = cgam  )

# text(x = da,
#      y = dg -0.1,
#     label = rownames(divpart.nat.perm$P) , col=cgam, cex =0.7 )
# 
# 
# text(x = da,
#      y = dg -0.1)
# ind <- locator(type="p")
text(x = ind$x,
     y = ind$y,
     label = rownames(divpart.nat.perm$P) , col=cgam, cex =0.7 )


abline(h=0, v=0, col="grey")
legend("topright", legend = c("native richness", "alien richness"), cex=0.85,
       pch= c(21,3), pt.bg=c("black","white"), bty="o")

# mtext(1,text = expression("SES"^alpha), line=2.5 )
# mtext(2,text = expression("SES"^gamma), line=2.5  , las=0)

 mtext(1,text = expression("SES ("*alpha*" richness)"), line=3 )
 mtext(2,text = expression("SES ("*gamma*" richness)"), line=2.5  , las=0)

# #### quantiles
# 
# dg<- (divpart.nat.perm$P$above.gamm -0.5)*2
# da<- (divpart.nat.perm$P$above.alpha - 0.5)*2
# dga<- (divpart.ali.perm$P$above.gamm -0.5)*2
# daa<- (divpart.ali.perm$P$above.alpha - 0.5)*2
# 
# x11()
# par( mar=c(4,4,2,2))
# plot(da, dg, xlim =c(-1,1),ylim =c(-1,1) , pch = 21, col=cgam ,bg = cgam, ann=F)
# points(daa, dga,  pch = 21, col=cgam )
# segments (da, dg, daa, dga,col=cgam)
# 
# text(x = da +  (c(0.00, 0.15, 0.00,  0, 0.0, 0.00, 0    ,0.0  , 0.00, 0.21, 0    )-0.2),
#      y = dg +    c(+0.1, -0.1, 0.00, 0 , 0.0, 0.00, -0.01, 0.0 , 0.00, -0.1, -0.03)
# )
# 
# ind <- locator()
# 
# text(x = ind$x,
#      y = ind$y,
#      label = rownames(divpart.nat.perm$P) , col=cgam, cex =0.7 )
# 
# abline(h=0, v=0, col="grey")
# legend("topright", legend = c("Native richness", "Alien Richness"), cex=0.9,
#        pch= c(21,21), pt.bg=c("black","white"), bty="n")
# 
# mtext(1,text = "Change in Alpha", line=2.5 )
# mtext(2,text = "Change in Gamma", line=2.5  )
