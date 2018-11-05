### Graphs and statistics for the paper on critical abundances

### Figures with NO covariabes:
#### Figure 2 NoCovar : trends in alpha richness effect size ####

M <- glmSRnat.overall.nocovar
sel = impsp.nocovar

par(mfrow = c(3,4), mar=c(0,0,2,1), oma=c(7,7,1,1))
ylim=c(0,150)
### Loop on selected species
for (i in 1:length(sel))  {
  # if (i ==4 | i == 7) plot.new()
  sp <- sel[i]
  es <- as.numeric(M$est[sp,])
  es <- exp(as.numeric(M$est[sp,]))*100
  n <- as.numeric(M$n.obs[sp,])[2:6]
  
  
  low <- as.numeric(M$CIlow[sp,])
  hi <- as.numeric(M$CIhi[sp,])
  low <- exp(as.numeric(M$CIlow[sp,]))*100
  hi <- exp(as.numeric(M$CIhi[sp,]))*100
  
  ## correct for small sample sizes
  small <- es
  small[which(n>=5)] <- NA
  es[which(n<5)] <- NA
  low[n<5] <- NA
  hi[n<5] <- NA
  
  # if ( !all(is.na(c(hi, low)))) lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
  # if ( all(is.na(c(hi, low)))) lims <-  max(c(abs(es)), na.rm=T) +0.02
  
  col <- rep("white",6)
  col[(1:6)>= as.numeric( M$impact.spread[sp,threshold])] <- "black"  # above threshold
  col[(2:6)[hi>100]] <- NA  # robust negative coef
  col[(2:6)[n<5]] <- NA  # sufficient data points
  
  # plot background
  plot(1:5,  rep(0,5), ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F, bty = "n")
  
  #dotted line
  # abline(h=0,lty="dotted")
  abline(h=exp(0)*100,lty="dotted")
  
  # draw small sample sizes
  # if (!all(is.na(small))) points(1:5, small, pch=21, cex = 0.6)
  
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
        font = 3, outer= F,adj=0.1, cex=0.7, line=0.2, las = 1)
  
  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  if (i %in% c(5:8)) text(y=-15, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
  
  # Y axis
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, labels = F, lwd = 0, lwd.ticks = 1)
  if ( i %in% c(1,5)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  
  # Y axis name
  if ( i %in% c(1,5)) {
    mtext(2, text=substitute("% SRnat"[a], list(a = "[rare]")), adj=0.5, cex=0.8, line=2, las = 0)
  }
  
  box(bty = "o", lwd = 1)
  
}

mtext(1, text=c("Abundance class"), adj=0.5, line=5, las = 1, outer=T)
mtext(2, text=expression(paste("Effect size on native ", alpha,"-richness")), adj=0.5, line=5, las = 0, outer=T)

#### Figure 3: trends in gamma richness effect size ####

par(mfrow = c(3,4), mar=c(0,0,2,1), oma=c(7,6,1,1))
ylim = c(-45,35)
xlim = c(0.5, 5.5)
for (i in 1 : length(sel)){
  
  # if (i ==4) plot.new()
  sp <- sel[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ (gamma.trend.nat$P.gamma[sp,] <=0.025) +1]
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
  
  # if (!all(is.na(small))) points(1:6, small, pch=21, cex = 0.6)
  
  # annotations
  
  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  # if (i %in% c(8:11)) text(y=-53.2, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
  if (i %in% c(5:8)) text(y=-53.2, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
  
  # Y axis
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, labels = F, lwd = 0, lwd.ticks = 1)
  # if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  if ( i %in% c(1,5)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  
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
        font = 3, outer= F,adj=0.1, cex=0.7, line=0.2, las = 1)
}

mtext(1, text=c("Abundance class"), adj=0.5, line=5, las = 1, outer=T)
mtext(2, text=expression(paste("Effect size on native ", gamma,"-richness")), adj=0.5, line=3, las = 0, outer=T)





#### __________________________________________DEFINE ELEMENTS__________________________________________ ####
abclasses= c("Rare" ,"Occasional",  "Frequent", "Common", "Abundant","Dominant")
threshold ="th.CI"
sel <- impsp
 # sel <- sel[c(1,9, 10,2:8)]

# loss in gamma native richness
# several different options of calculation but we choose the first one because it is a statistic on the difference (below - above) itself (like a z-statistic for two sample test) rather than a statistics on the mean above (= one sample test).
dg <- table.div.part$deltagamma - table.div.part$delta.null.permute.all # observed difference in means minus the expected difference in means (= sort of z statistics for comparing means of above and below gamma) ?
dg2 <- table.div.part$delta.z.permute.all # observed difference in means 
dg0<- table.div.part$gamma.loss  ### Difference between observed and expected mean above crit.
cor.test(dg, dg0)
plot(-dg, dg0, xlab = "delta.gamma.c", type = "n", ylab = "gamma.above.c", ylim = c(0,60), xlim = c(0,60))
abline(0,1)
text(-dg, dg0, labels = rownames(table.div.part), cex = 0.7)

dg1<- table.div.part$deltagamma - table.div.part$delta.null.beta
dg1<- table.div.part$delta.z.beta


# loss in alpha native richness
da<- table.div.part$aRc -  table.div.part$aRo
da1<- table.div.part$aRc -  table.div.part$aRnull
da2 <- table.div.part$deltaalpha.z.permute.all

plot(da, da2)
abline(0,1)

# loss in alpha native richness
dbeta<- table.div.part$z.beta.diff


#spatial distributions:
spread <- table.div.part$n.plot.impact
prevalence <- table.div.part$prevalence
dominance <- table.div.part$n.plot.dominant
prop.impact <- spread/prevalence

# critical abundance level
critical.abun <- table.div.part$th.CI

# library(FactoMineR)
# tmp <- PCA(data.frame(cbind(dg0, dg1,dg, da,da1, dbeta), row.names= impsp ) )
# plot(tmp, choix = "ind")
# plot(tmp, choix = "var")
# plot(tmp, choix = "var", axes = c(1,3))

### Looking at GLM results across species and covariables
names(glmSRnat.overall)[6] <- "covar.tab"
# Distribution of GLM coefs
boxplot(cbind(Elevation = glmSRnat.overall$covar.tab$DEM_10$coef.glm,
              Slope = glmSRnat.overall$covar.tab$SLOPE$coef.glm, 
              Northness = glmSRnat.overall$covar.tab$Northern$coef.glm,
              SRali = glmSRnat.overall$covar.tab$SRali$coef.glm, 
              FocalSpecies = rowMeans(glmSRnat.overall$est)),
        outline = FALSE, ylab = "GLM coefficients")
abline(h = 0, col = "lightgrey")

# Frequency of significance for each factor:

barplot(cbind(Elevation = table(glmSRnat.overall$covar.tab$DEM_10$P.coef<0.05),
              Slope =  table(glmSRnat.overall$covar.tab$SLOPE$P.coef<0.05), 
              Northness =  table(glmSRnat.overall$covar.tab$Northern$P.coef<0.05),
              SRali =  table(glmSRnat.overall$covar.tab$SRali$P.coef<0.05) ),
        ylab = "GLM signif")





####### Calculate correlations within table 1   #########
library(corrplot)
mat <- table2[, c(1:4, 9, 13, 15)] 
# # log option:
# mat$Presence <- log(mat$Presence)
# mat$above.Acrit <- log (mat$above.Acrit)
# mat$dominance <- log (mat$dominance+1)

cor.mtest <- function(mat, conf.level = 0.95, method = "pearson"){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- coef.mat <- df.mat <- matrix(NA, n, n)
  diag(p.mat) = 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], method = method)
      
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      coef.mat[i,j] <- coef.mat[j,i] <- tmp$estimate
      if (!is.null(tmp$parameter)) df.mat[i,j] <- df.mat[j,i] <- tmp$parameter
    }
  }
  return(list(coef = coef.mat, df = df.mat, p = p.mat ))
}
res1 <- cor.mtest(mat)
res2 <- cor.mtest(mat, method = "spearman")

## plot
par(mfrow = c(1,2), mar = c(1,5,8,2))
corrplot(cor(mat), p.mat = res1[[3]], sig.level=0.05, insig = "blank")
title(main = "Pearson")
corrplot(cor(mat, method = "spearman"), p.mat = res2[[3]], sig.level=0.05, insig = "blank")
title(main = "Spearman")

# export the correlation tables in csv

pearson.mat <- NULL
for (i in 1: length(names(mat))) {
  tmp <- rbind(coef = round(res1[[1]][i,], 4),
               P.val = round(res1[[3]][i,], 4),
               df = res1[[2]][i,])
  pearson.mat = rbind(pearson.mat, tmp)
}
colnames(pearson.mat) = names(mat)

spearman.mat <- NULL
for (i in 1: length(names(mat))) {
  tmp <- rbind(coef = round(res2[[1]][i,],4),
               P.val = round(res2[[3]][i,], 4),
               df = res2[[2]][i,])
  spearman.mat = rbind(spearman.mat, tmp)
}
colnames(spearman.mat) = names(mat)

write.csv(pearson.mat, "Pearson correlations.csv")
write.csv(spearman.mat, "Spearman correlations.csv")


##colour code For figures
# cgam <- c("grey50", "black")[ (table.div.part$GRP.permute.all<=0.05) + 1]
cgam <- "black"








## Fig 5 alternative #########
#loss in alpha vs. loss in gamma (lg transformed pearson corr)   
par(mfrow = c(1,1), mar=c(3,1,2,1), oma=c(1,2,0,0), las = 1, xpd = TRUE)
plot(da, dg ,pch = 20, col=cgam , ann=F,type ="n", axes =F, log = '', xlim = c(-4,0))
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x = da,y = dg,rownames(table.div.part) , col="grey", cex =0.6, pos = 4)
points(x = da,y = dg, pch = 20)
box(bty="l")
f <- cor.test(da ,dg, method = "pearson")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=F))),adj = 1, cex=0.8, font = 3, line= 0.5)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1, line= 0.5)
mtext(1,text = expression(Delta*alpha*"-richness"), line=2 )
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=2 , las=0)




#### Figure 4 alternative: Maps + correlations in one figure  #######
db = databp[databp$PlotName %in% unimprovedgrasslands, ]
targets <- sel

par(
  mfrow = c(4,3),
  mar = c(0, 0, 1, 0),
  oma = c(2, 2, 0, 0)
)



for (i in 1:9) {
  if (i ==9) {
    plot(region, border = NA)
    legend(
      'center',
      legend = c("< critical abundance", "> critical abundance", "dominant"),
      fill = c("tan", "sienna2", "sienna4"),
      xpd = NA,
      bty = "n",
      border = c("tan", "sienna2", "sienna4")
    )
    addnortharrow(pos = "topright", padin = c(0.40, 0), scale = 0.5,
                  lwd = 1, border = "black", cols = c("white", "black"),
                  text.col = "black")
    map.scale(x = extent(region)[1]+0.1, y = extent(region)[3] + 0.08, ratio = F, relwidth = 0.40)
  }
  
  if (i <9) {
    sp=targets[i]
    plot.alien = as.character(db[which(db$abun %in% c(1, 2, 3, 4, 5, 6) &
                                         db$vegtype == "G" & db$SpeciesCode == sp), "PlotName"])
    th = glmSRnat.overall$impact.spread[sp, "th.CI"]
    es <- as.numeric(glmSRnat.overall$est[sp, ])
    n <- as.numeric(glmSRnat.overall$n.obs[sp, ])[2:6]
    es[which(n < 5)] <- NA
    es = which(!is.na(es)) + 1
    m <- max(es, na.rm = T)
    max.impact = as.character(db[which(db$abun >= 6 &
                                         db$vegtype == "G" & db$SpeciesCode == sp), "PlotName"])
    th.impact = as.character(db[which(db$abun >= th &
                                        db$vegtype == "G" & db$SpeciesCode == sp), "PlotName"])
    
    
    plot(region, col = "grey50", border = NA)
    polygon(extent(region)[c(1,2,2,1)], extent(region)[c(3,3,4,4)], col ="lightblue", border = NA)
    plot(region, add= T, col = "grey70", border = NA)
    plot(study_area, add= T, col = "grey90", border = NA)
    plot(region, add= T, border = "grey70")
    
    plot(
      envplot[plot.alien,],
      pch = 22,
      cex = 0.5,
      col = "tan",
      bg = "tan",
      add = T
    )
    plot(
      envplot[th.impact,],
      pch = 22,
      cex = 0.5,
      col = "sienna2",
      bg = "sienna2",
      add = T
    )
    plot(
      envplot[max.impact,],
      pch = 22,
      cex = 0.5,
      col = "sienna4",
      bg = "sienna4",
      add = T
    )
    
    mtext(
      3,
      text = paste(letters[i],") ",sub("_", " ", species[sp, "tip"]), sep = ""),
      adj = 0.1,
      line = -0.5,
      font = 3,
      cex = 0.6
    )
  }
  
}

## underneath: add the correlations with gamma loss
dg.lab = dg
dg.lab [ rownames(table.div.part) == "ANTODO"] <- -23
dg.lab [ rownames(table.div.part) == "DACGLO"] <- -26
# 
par( mar=c(2,1,3,1),  las = 1)


### Prevalence
print(f <- cor.test(log(prevalence),dg, method = "spearman"))
# print(f <- cor.test(prevalence[-c(1,8)],dg[-c(1,8)]))
plot (prevalence,dg , pch = 20, ann=F, log = "x", xlim =c(150,850), type= "p", axes=F)
text(x = prevalence, y = dg.lab, rownames(table.div.part) , col="tan", cex =0.7, pos =4)
points(x = prevalence,y = dg, pch = 20, col="tan")
box(bty="l")
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
box(bty="l")
mtext(3, text = substitute(rho == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.7,line = 0.5)
mtext(1,text = "Presence", line=1.5, cex = 0.7, font =3)
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=1.5, las=0, cex = 0.8)
mtext(3, text = 'i)', adj = 0,
      line = 0.5,
      font = 3,
      cex = 0.6)



### spread above critical abundance
print(f <- cor.test(log(spread), dg, method = "spearman"))
plot (spread, dg , pch = 20, ann=F, log = "x", xlim =c(5, 2000), type ="p", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=F, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x = spread, y = dg, rownames(table.div.part) , col="sienna2", cex =0.7, pos = 4)
points(x = spread,y = dg, pch = 20, col="sienna2")
box(bty="l")
mtext(3, text = substitute(italic(rho == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.7,line = 0.5)
mtext(1,text = "> critical abundance", line=1.5 , cex = 0.7, padj = 0, font =3)
mtext(3, text = 'j)', adj = 0,
      line = 0.5,
      font = 3,
      cex = 0.6)



### dominance
print(f <- cor.test(dominance,dg, method= "spearman"))
plot (dominance+1, dg , pch = 20, ann=F, log = "x", type= "p",axes =F, xlim = c(1,100))
text(x = dominance+1, y = dg.lab, rownames(table.div.part) , col="sienna4", cex =0.7, pos = 4)
points(x = dominance+1,y = dg, pch = 20, col="sienna4")
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=F, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
box(bty="l")
mtext(3, text = substitute(rho == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.7,line = 0.5)
mtext(1,text = "Dominance", line=1.5, cex = 0.7, font =3)
mtext(3, text = 'k)', adj = 0,
      line = 0.5,
      font = 3,
      cex = 0.6)

mtext(1,text = "Number of plots", line=1, cex = 0.8, outer =TRUE)









###___________________________=___________ EXTRA FIGURES AND STATS ______________________________________####

## Extra figure for reviewers showing no longer significant 3 species :

par(mfrow = c(1,3), mar=c(0,0,2,1), oma=c(7,7,1,1))

M <- glmSRnat.overall
nonsigsp <- impsp.nocovar[which(! impsp.nocovar %in% impsp)]
threshold = "th.CI"
ylim=c(0,200)
### Loop on selected species
for (i in 1:length(nonsigsp))  {
  # if (i ==4 | i == 7) plot.new()
  sp <- nonsigsp[i]
  es <- as.numeric(M$est[sp,])
  es <- exp(as.numeric(M$est[sp,]))*100
  n <- as.numeric(M$n.obs[sp,])[2:6]
  
  low <- as.numeric(M$CIlow[sp,])
  hi <- as.numeric(M$CIhi[sp,])
  low <- exp(as.numeric(M$CIlow[sp,]))*100
  hi <- exp(as.numeric(M$CIhi[sp,]))*100
  
  ## correct for small sample sizes
  small <- es
  small[which(n>=5)] <- NA
  es[which(n<5)] <- NA
  low[n<5] <- NA
  hi[n<5] <- NA
  
  # if ( !all(is.na(c(hi, low)))) lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
  # if ( all(is.na(c(hi, low)))) lims <-  max(c(abs(es)), na.rm=T) +0.02
  
  col <- rep("white",6)
  col[(1:6)>= as.numeric( M$impact.spread[sp,threshold])] <- "black"  # above threshold
  col[(2:6)[hi>100]] <- NA  # robust negative coef
  col[(2:6)[n<5]] <- NA  # sufficient data points
  
  # plot background
  plot(1:5,  rep(0,5), ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F, bty = "n")
  
  #dotted line
  # abline(h=0,lty="dotted")
  abline(h=exp(0)*100,lty="dotted")
  
  # draw small sample sizes
  # if (!all(is.na(small))) points(1:5, small, pch=21, cex = 0.6)
  
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
        font = 3, outer= F,adj=0.1, cex=0.7, line=0.2, las = 1)
  
  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  text(y=-15, x = 1:5, labels= abclasses[2:6],  cex=1, srt=60, adj=1, xpd = NA)
  
  # Y axis
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, labels = F, lwd = 0, lwd.ticks = 1)
  if ( i %in% c(1)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  
  # Y axis name
  if ( i %in% c(1)) {
    mtext(2, text=substitute("% SRnat"[a], list(a = "[rare]")), adj=0.5, cex=0.8, line=2, las = 0)
  }
  
  box(bty = "o", lwd = 1)
  
}

mtext(1, text=c("Abundance class"), adj=0.5, line=5, las = 1, outer=T)
mtext(2, text=expression(paste("Effect size on native ", alpha,"-richness")), adj=0.5, line=5, las = 0, outer=T)



## figure S2:  trends in alpha richness effect size without covar ####

par(mfrow = c(3,4), mar=c(0,0,2,1), oma=c(7,7,1,1))

M <- glmSRnat.overall.nocovar
ylim=c(0,250)
### Loop on selected species
for (i in 1:length(impsp.nocovar))  {
  # if (i ==4 | i == 7) plot.new()
  sp <- impsp.nocovar[i]
  es <- as.numeric(M$est[sp,])
  es <- exp(as.numeric(M$est[sp,]))*100
  n <- as.numeric(M$n.obs[sp,])[2:6]
  
  low <- as.numeric(M$CIlow[sp,])
  hi <- as.numeric(M$CIhi[sp,])
  low <- exp(as.numeric(M$CIlow[sp,]))*100
  hi <- exp(as.numeric(M$CIhi[sp,]))*100
  
  ## correct for small sample sizes
  small <- es
  small[which(n>=5)] <- NA
  es[which(n<5)] <- NA
  low[n<5] <- NA
  hi[n<5] <- NA
  
  # if ( !all(is.na(c(hi, low)))) lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
  # if ( all(is.na(c(hi, low)))) lims <-  max(c(abs(es)), na.rm=T) +0.02
  
  col <- rep("white",6)
  col[(1:6)>= as.numeric( M$impact.spread[sp,threshold])] <- "black"  # above threshold
  col[(2:6)[hi>100]] <- NA  # robust negative coef
  col[(2:6)[n<5]] <- NA  # sufficient data points
  
  # plot background
  plot(1:5,  rep(0,5), ylim=ylim, xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F, bty = "n")
  
  #dotted line
  # abline(h=0,lty="dotted")
  abline(h=exp(0)*100,lty="dotted")
  
  # draw small sample sizes
  # if (!all(is.na(small))) points(1:5, small, pch=21, cex = 0.6)
  
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
        font = 3, outer= F,adj=0.1, cex=0.7, line=0.2, las = 1)
  
  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  if (i %in% c(8:11)) text(y=-15, x = 1:5, labels= abclasses[2:6],  cex=1, srt=90, adj=1, xpd = NA)
  
  # Y axis
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, labels = F, lwd = 0, lwd.ticks = 1)
  if ( i %in% c(1,5,9)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  
  # Y axis name
  if ( i %in% c(1,5,9)) {
    mtext(2, text=substitute("% SRnat"[a], list(a = "[rare]")), adj=0.5, cex=0.8, line=2, las = 0)
  }
  
  box(bty = "o", lwd = 1)
  
}

mtext(1, text=c("Abundance class"), adj=0.5, line=5, las = 1, outer=T)
mtext(2, text=expression(paste("Effect size on native ", alpha,"-richness")), adj=0.5, line=5, las = 0, outer=T)




#### Figure 6 alternative: Maps of presence and impact of two selected species ####

db = databp[databp$PlotName %in% unimprovedgrasslands, ]
targets <- c('ACHMIL', 'CYNCRI')

par(
  mfrow = c(1, 3),
  mar = c(0, 0, 0, 0),
  oma = c(2, 0, 0, 0)
)

for (i in targets) {
  plot.alien = as.character(db[which(db$abun %in% c(1, 2, 3, 4, 5, 6) &
                                       db$vegtype == "G" & db$SpeciesCode == i), "PlotName"])
  th = glmSRnat.overall$impact.spread[i, "th.CI"]
  es <- as.numeric(glmSRnat.overall$est[i, ])
  n <- as.numeric(glmSRnat.overall$n.obs[i, ])[2:6]
  es[which(n < 5)] <- NA
  es = which(!is.na(es)) + 1
  m <- max(es, na.rm = T)
  max.impact = as.character(db[which(db$abun >= 6 &
                                       db$vegtype == "G" & db$SpeciesCode == i), "PlotName"])
  th.impact = as.character(db[which(db$abun >= th &
                                      db$vegtype == "G" & db$SpeciesCode == i), "PlotName"])
  
  plot(study_area, col = "white", border = "grey60")
  plot(
    envplot[plot.alien,],
    pch = 22,
    cex = 0.5,
    col = "tan",
    bg = "tan",
    add = T
  )
  plot(
    envplot[th.impact,],
    pch = 22,
    cex = 0.5,
    col = "sienna2",
    bg = "sienna2",
    add = T
  )
  plot(
    envplot[max.impact,],
    pch = 22,
    cex = 0.5,
    col = "sienna4",
    bg = "sienna4",
    add = T
  )
  
  mtext(
    3,
    text = sub("_", " ", species[i, "tip"]),
    adj = 0.1,
    line = -1.5,
    font = 3,
    cex = 0.7
  )
}

plot.new()
legend(
  'center',
  legend = c("< critical abundance", ">critical abundance", "Dominant"),
  fill = c("tan", "sienna2", "sienna4"),
  xpd = NA,
  bty = "n",
  border = c("tan", "sienna2", "sienna4")
) 



## Three species
db = databp[databp$PlotName %in% unimprovedgrasslands, ]
targets <- c('ACHMIL', 'CYNCRI', "ANTODO")

par(
  mfrow = c(1, 3),
  mar = c(0, 0, 0, 0),
  oma = c(2, 0, 0, 0)
)

for (i in targets) {
  plot.alien = as.character(db[which(db$abun %in% c(1, 2, 3, 4, 5, 6) &
                                       db$vegtype == "G" & db$SpeciesCode == i), "PlotName"])
  th = glmSRnat.overall$impact.spread[i, "th.CI"]
  es <- as.numeric(glmSRnat.overall$est[i, ])
  n <- as.numeric(glmSRnat.overall$n.obs[i, ])[2:6]
  es[which(n < 5)] <- NA
  es = which(!is.na(es)) + 1
  m <- max(es, na.rm = T)
  max.impact = as.character(db[which(db$abun >= 6 &
                                       db$vegtype == "G" & db$SpeciesCode == i), "PlotName"])
  th.impact = as.character(db[which(db$abun >= th &
                                      db$vegtype == "G" & db$SpeciesCode == i), "PlotName"])
  
  plot(study_area, col = "white", border = "grey60")
  plot(
    envplot[plot.alien,],
    pch = 22,
    cex = 0.5,
    col = "tan",
    bg = "tan",
    add = T
  )
  plot(
    envplot[th.impact,],
    pch = 22,
    cex = 0.5,
    col = "sienna2",
    bg = "sienna2",
    add = T
  )
  plot(
    envplot[max.impact,],
    pch = 22,
    cex = 0.5,
    col = "sienna4",
    bg = "sienna4",
    add = T
  )
  
  mtext(
    3,
    text = sub("_", " ", species[i, "tip"]),
    adj = 0.1,
    line = -1.5,
    font = 3,
    cex = 0.7
  )
} 

# #### OLD Figure 6: Maps of presence and impact ####
# 
# db = databp[databp$PlotName %in% unimprovedgrasslands, ]
# targets <- sel
# 
# par(
#   mfrow = c(4, 3),
#   mar = c(0, 0, 1, 0),
#   oma = c(2, 0, 0, 0),
#   xpd = T
# )
# 
# for (i in 1:12) {
#   if (i ==12) {
#     plot(region, border = NA)
#     legend(
#       'center',
#       legend = c("< critical abundance", "> critical abundance", "dominant"),
#       fill = c("tan", "sienna2", "sienna4"),
#       xpd = NA,
#       bty = "n",
#       border = c("tan", "sienna2", "sienna4")
#     ) 
#     
#       map.scale(x = extent(region)[1]+0.2, y = extent(region)[3] + 0.05, ratio = F) 
#   }
#   if (i <12) {
#   sp=targets[i]
#   plot.alien = as.character(db[which(db$abun %in% c(1, 2, 3, 4, 5, 6) &
#                                        db$vegtype == "G" & db$SpeciesCode == sp), "PlotName"])
#   th = glmSRnat.overall$impact.spread[sp, "th.CI"]
#   es <- as.numeric(glmSRnat.overall$est[sp, ])
#   n <- as.numeric(glmSRnat.overall$n.obs[sp, ])[2:6]
#   es[which(n < 5)] <- NA
#   es = which(!is.na(es)) + 1
#   m <- max(es, na.rm = T)
#   max.impact = as.character(db[which(db$abun >= 6 &
#                                        db$vegtype == "G" & db$SpeciesCode == sp), "PlotName"])
#   th.impact = as.character(db[which(db$abun >= th &
#                                       db$vegtype == "G" & db$SpeciesCode == sp), "PlotName"])
#   
#   plot(study_area, col = "white", border = "grey60")
#   plot(
#     envplot[plot.alien,],
#     pch = 22,
#     cex = 0.5,
#     col = "tan",
#     bg = "tan",
#     add = T
#   )
#   plot(
#     envplot[th.impact,],
#     pch = 22,
#     cex = 0.5,
#     col = "sienna2",
#     bg = "sienna2",
#     add = T
#   )
#   plot(
#     envplot[max.impact,],
#     pch = 22,
#     cex = 0.5,
#     col = "sienna4",
#     bg = "sienna4",
#     add = T
#   )
#   
#   mtext(
#     3,
#     text = paste(letters[i],") ",sub("_", " ", species[sp, "tip"]), sep = ""),
#     adj = 0.1,
#     line = 0,
#     font = 3,
#     cex = 0.6
#   )
#   }
# }
# 


#### Compare 3 different metrics of spatial spread vs. gamma loss #####

#GLM model for explaining deltagamma => not better
summary(f1 <- glm((-dg) ~ prevalence, family = quasipoisson(link = "log") ))
summary(f2 <- glm((-dg) ~ spread, family = quasipoisson(link = "log") ))
summary(f3 <- glm((-dg) ~ dominance, family = quasipoisson(link = "log") ))
summary( glm( -dg ~ da , family = quasipoisson(link = "log") ))
summary( glm( -dg ~ da*spread , family = quasipoisson(link = "log") ))
summary( glm( -dg ~ da , family = quasipoisson(link = "log") ))


#GLM model for explaining gamma.corrected
print(f <- cor.test(prevalence,dg1))
plot(prevalence,dg1)
print(f <- cor.test(spread,dg1))
plot(spread, dg1)
print(f <- cor.test(dominance,dg1))


par(mfrow = c(1,3), mar=c(3,3,2,1), oma=c(1,2,0,0), las = 1)
### Prevalence
print(f <- cor.test(prevalence,dg, method = "spearman"))
# print(f <- cor.test(prevalence[-c(1,8)],dg[-c(1,8)]))
plot (prevalence,dg , pch = 20, ann=F, log = "", xlim =c(100, 600), type= "p", axes=F)
text(prevalence,dg, label = impsp)
# text(x = prevalence, y = dg,label = rownames(table.div.part),  col=cgam ,cex=0.7)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
box(bty="l")
mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.8)
mtext(1,text = "Number of plots", line=2)
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=2 , las=0)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1)

### spread above critical abundance
print(f <- cor.test(spread ,log(-dg)))
plot (spread, dg , pch = 20, ann=F, log = "", xlim =c(2.5, 600), type ="p", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=F, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
# text(x = spread, y = dg,label = rownames(table.div.part),col=cgam, cex=0.7)

box(bty="l")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8)
mtext(1,text = "Number of plots > critical abundance", line=2 )
mtext(3, text = 'b)',adj = 0, cex=0.8, font = 1)


### proportion of spread above critical abundance
print(f <- cor.test(prop.impact ,dg))
plot (prop.impact, dg , pch = 20, ann=F, log = "",  type ="p", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=F, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
# text(x = spread, y = dg,label = rownames(table.div.part),col=cgam, cex=0.7)

box(bty="l")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8)
mtext(1,text = "% of plots > critical abundance", line=2 )
mtext(3, text = 'c)',adj = 0, cex=0.8, font = 1)

### dominance
print(f <- cor.test(dominance,dg))
plot (dominance, dg , pch = 20, ann=F, log = "", type= "p",axes =F)
# text(x = dominance, y = dg,label = rownames(table.div.part),  col=cgam ,cex=0.7)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=F, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
box(bty="l")
mtext(3, text = substitute(r == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.8)
mtext(1,text = "Number of dominated plots", line=2)
mtext(3, text = 'c)',adj = 0, cex=0.8, font = 1)


#### Loss in beta vs. loss in gamma ####
par(mfrow = c(1,2), mar=c(3,1,2,1), oma=c(1,2,0,0), las = 1, xpd = TRUE)
plot(db, dg ,pch = 20,ann=F, type ="n", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x = db,y = dg,rownames(table.div.part) , col="grey", cex =0.6, pos = 4)
points(x = db,y = dg, pch = 20)
box(bty="l")
 f <- cor.test(db ,dg, method = "pearson")
 mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=F))),adj = 1, cex=0.8, font = 3, line= 0.5)
# 
# f <- cor.test(db[-10] ,dg[-10], method = "pearson")
# mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=F))),adj = 1, cex=0.8, font = 3, line= 0.5)

mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1, line= 0.5)
mtext(1,text = expression(Delta*beta*"-diversity"), line=2 )
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=2 , las=0)

plot(spread, db ,pch = 20, col=cgam , ann=F, type ="n", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x = spread,y = db,rownames(table.div.part) , col="grey", cex =0.6, pos = 4)
points(x = spread,y = db, pch = 20)
box(bty="l")
f <- cor.test(db ,spread, method = "spearman")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=F))),adj = 1, cex=0.8, font = 3, line= 0.5)

# f <- cor.test(db[-10] ,spread[-10], method = "pearson")
# mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=F))),adj = 1, cex=0.8, font = 3, line= 0.5)

mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1, line= 0.5)
mtext(2,text = expression(Delta*beta*"-diversity"), line=2, las = 0 )
mtext(1,text = expression("Distribution of potential impact"), line=2 )


#### Loss in alpha corrected vs. loss in gamma corrected by NM ####
par(mfrow = c(1,2), mar=c(3,1,2,1), oma=c(1,2,0,0), las = 1, xpd = TRUE)
plot(da2, dg ,pch = 20,ann=F, type ="n", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x = da2,y = dg,rownames(table.div.part) , col="grey", cex =0.6, pos = 4)
points(x = da2,y = dg, pch = 20)
box(bty="l")
f <- cor.test(da2 ,dg, method = "pearson")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=F))),adj = 1, cex=0.8, font = 3, line= 0.5)
# 
# f <- cor.test(db[-10] ,dg[-10], method = "pearson")
# mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=F))),adj = 1, cex=0.8, font = 3, line= 0.5)

mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1, line= 0.5)
mtext(1,text = expression(Delta*beta*"-diversity"), line=2 )
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=2 , las=0)


#### Spatial distribution of species #####
par(mfrow = c(1,2), mar=c(3,3,2,1), oma=c(1,2,0,0), las = 1)

print(f <- cor.test(prevalence,(spread), method = "pearson"))
plot (prevalence,spread , pch = 20, ann=F, type= "p", axes=F, log = "")
# text(x = prevalence, y = dg,label = rownames(table.div.part),  col=cgam ,cex=0.7)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
box(bty="l")
mtext(3, text = substitute(rho == est *" "*p, list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),
      adj = 1, cex=0.8)
mtext(1,text = "Prevalence", line=2)
mtext(2,text = "> critical abundance", line=2 , las=0)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1)

### spread above critical abundance
print(f <- cor.test(spread , dominance, method = "spearman"))
plot (dominance, spread , pch = 20, col="black" , ann=F, log = "", type ="p", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=F, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
# text(x = spread, y = dg,label = rownames(table.div.part),col=cgam, cex=0.7)
box(bty="l")
mtext(3, text = substitute(italic(rho == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8)
mtext(1,text = "Plots (species = dominant)", line=2 )
mtext(3, text = 'b)',adj = 0, cex=0.8, font = 1)


#### Critical abundance levels vs. magnitude of impact #####
par(mfrow = c(1,2), mar=c(3,2,2,1), oma=c(1,2,0,0), las = 1)

#loss in alpha vs. critical abun

plot(critical.abun, da2 ,pch = 20, col=cgam , ann=F,type ="n", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x =critical.abun ,y = da2,label = rownames(table.div.part) , col=cgam, cex =0.7)
f <- cor.test(da2 ,critical.abun)
box(bty="o")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8, font = 3)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1)
mtext(2,text = expression(Delta*alpha*"-richness"), line=2, las = 0)
mtext(1,text = "Critical abundance", line=2 )


#loss in gamma vs. critical abun
plot(critical.abun, dg ,pch = 20, col=cgam , ann=F,type ="n", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x =critical.abun ,y = dg,label = rownames(table.div.part) , col=cgam, cex =0.7)
f <- cor.test(dg ,critical.abun)
box(bty="o")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8, font = 3)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1)
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line= 2, las = 0)
mtext(1,text = "Critical abundance", line=2 )



#### Critical abundance levels vs. spatial spread #####
par(mfrow = c(1,3), mar=c(3,2,2,1), oma=c(1,2,0,0), las = 1)

# Prevalence vs. critical abun
plot(critical.abun, prevalence ,pch = 20, col=cgam , ann=F,type ="n", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x =critical.abun ,y = prevalence,label = rownames(table.div.part) , col=cgam, cex =0.7)
f <- cor.test(prevalence ,critical.abun)
box(bty="o")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8, font = 3)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1)
mtext(2,text = "prevalence", line=2, las = 0)
mtext(1,text = "Critical abundance", line=2 )


# Spread of potential impact vs. critical abun
plot(critical.abun, spread ,pch = 20, col=cgam , ann=F,type ="n", log = "y", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x =critical.abun ,y = spread,label = rownames(table.div.part) , col=cgam, cex =0.7)
f <- cor.test(log(spread) ,critical.abun)
box(bty="o")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8, font = 3)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1)
mtext(2,text = "spread", line=2, las = 0)
mtext(1,text = "Critical abundance", line=2 )


# dominance vs. critical abun
plot(critical.abun, dominance ,pch = 20, col=cgam , ann=F,type ="n", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x =critical.abun ,y = dominance,label = rownames(table.div.part) , col=cgam, cex =0.7)
f <- cor.test(dominance ,critical.abun)
box(bty="o")
mtext(3, text = substitute(italic(r == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8, font = 3)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1)
mtext(2,text = "dominance", line=2, las = 0)
mtext(1,text = "Critical abundance", line=2 )


#######  Counting number of alien sp per plot  ########
library(ape)
db<- databp[databp$PlotName %in% unimprovedgrasslands,]

focal.com <- table(db[db$SpeciesCode %in% sel,c('SpeciesCode','PlotName')])
focal.com <- focal.com[,colSums(focal.com) >0]

plot(hclust(dist(focal.com)))
heatmap(focal.com)

nb.focal <- table(db[db$SpeciesCode %in% sel,'PlotName'])
focal.com <- comm[unimprovedgrasslands, sel]

hist(rowSums(focal.com>1))


focal.frequency <- as.data.frame(rowSums(ceiling(focal.com/10)))
colnames(focal.frequency) <- "presence"
focal.frequency$above.critical <- colSums(apply(focal.com, 1, function(x) {
  x >= glmSRnat.overall$impact.spread[sel, "th.CI"]
}))
focal.frequency$dom <- colSums(apply(focal.com, 1, function(x) {
  x >= 5
}))

focal.frequency$SR <- envplot$SR [rownames(focal.frequency)]
focal.frequency$SRnat <- envplot$SRnat [rownames(focal.frequency)]
focal.frequency$SRali <- envplot$SRali [rownames(focal.frequency)]

par(mfrow=c(2,2), mar = c(3,3,1,1), oma =c(4,3,1,1), xpd = NA, las = 1)
plot(SRnat ~ jitter(presence,1), focal.frequency, col = "#7F7F7F80", xlab="", ylab =expression(paste("Native ", alpha,"-richness")))
f <- cor.test(focal.frequency$SRnat ,as.numeric(focal.frequency$presence), method = "spearman")
mtext(side= 3, text =substitute(rho*" = "*r*p, list(r=round(f$estimate,2), p = p2star(f$p.value))), adj = 1, cex = 0.8)

plot(SRnat ~ jitter(above.critical, 1), focal.frequency, col = "#7F7F7F80",ylab = "", xlab="")
f <- cor.test(focal.frequency$SRnat ,as.numeric(focal.frequency$above.critical), method = "spearman")
mtext(side= 3, text =substitute(rho*" = "*r*p, list(r=round(f$estimate,2), p = p2star(f$p.value))), adj = 1, cex = 0.8)

plot(SRali ~ jitter(presence,1), focal.frequency, col = "#7F7F7F80", xlab = "",  ylab =expression(paste("Alien ", alpha,"-richness")))
mtext(1, text = "Number of co-occurring significant alien species", line= 3.5, cex = 0.8)

f <- cor.test(focal.frequency$SRali ,as.numeric(focal.frequency$presence), method = "spearman")
mtext(side= 3, text =substitute(rho*" = "*r*p, list(r=round(f$estimate,2), p = p2star(f$p.value))), adj = 1, cex = 0.8)

plot(SRali  ~ jitter(above.critical, 1), focal.frequency, col = "#7F7F7F80",  ann = F)
mtext(1, text = "Number of significant alien species\n co-occurring above critical abundance", line= 3.5, cex = 0.8)

f <- cor.test(focal.frequency$SRali ,as.numeric(focal.frequency$above.critical), method = "spearman")
mtext(side= 3, text =substitute(rho*" = "*r*p, list(r=round(f$estimate,2), p = p2star(f$p.value))), adj = 1, cex = 0.8)

par(mfrow=c(1,2))
plot(SR ~  jitter(presence,1), focal.frequency)
plot(SR ~ jitter(above.critical, 1), focal.frequency)


#### Figure 3: trends in BETA PROP ####
quartz()
par(mfrow = c(4,6), mar = c(0,0,2,0))
sapply(aliens[aliens%in% rownames(beta.trend.nat$obs)], function(i) {
  n <- na.omit(glmSRnat.overall$n.obs[i,])
  x <- beta.trend.nat$obs[i,]
  x[n < 5] <- NA
  plot(1:length(x),x, type = "b", xlim = c(1,7), ylim = c(0,1), axes = F)
  box(bty = "l")
  mtext(1, at =  1:length(x),line = -1, text = na.omit(glmSRnat.overall$n.obs[i,]), cex = 0.6) 
  mtext(3, text = i, cex = 0.4) 
  })

par(mfrow = c(4,6), mar = c(0,0,2,0))
sapply(natives[natives %in% rownames(beta.trend.nat$obs)], function(i) {
  n <- na.omit(glmSRnat.overall$n.obs[i,])
  x <- beta.trend.nat$obs[i,]
  x[n < 5] <- NA
  plot(1:length(x),x, type = "b", xlim = c(1,7), ylim = c(0,1), axes = F)
  box(bty = "l")
  mtext(1, at =  1:length(x),line = -1, text = na.omit(glmSRnat.overall$n.obs[i,]), cex = 0.6) 
  mtext(3, text = i, cex = 0.4) 
})
