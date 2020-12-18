### Graphs and statistics for the paper on critical abundances


#### __________________________________________DEFINE ELEMENTS__________________________________________ ####
abclasses= c("Rare" ,"Occasional",  "Frequent", "Common", "Abundant","Dominant")
threshold ="th.CI"
sel <- impsp
 # sel <- sel[c(1,9, 10,2:8)]
table.div.part <- table.div.part[sel,]
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

# Recalculate loss in alpha at th
da <- sapply(sel, FUN = function(i) {
  th <- table.div.part[i,"th.CI"]
  return(alpha.above.trend$alpha.above[i,th] -alpha.above.trend$alpha.below[i,th]) }
)
# max loss in alpha
da.max <- sapply(sel, FUN = function(i) {
  m<- max(which(!is.na(glmSRnat.overall$est[i,]) & (glmSRnat.overall$n.obs[i,2:6] >= min.occur)), na.rm = T)+1
  return(alpha.above.trend$alpha.above[i,m] -alpha.above.trend$alpha.below[i,m]) }
  )

ab.max <- sapply(sel, FUN = function(i) {
  m<- max(which(!is.na(glmSRnat.overall$est[i,]) & (glmSRnat.overall$n.obs[i,2:6] >= min.occur)), na.rm = T)+1
  return(m) }
)


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

par(mfrow= c(1,2))
# Distribution of GLM coefs
boxplot(cbind(Year = exp(glmSRnat.overall$covar.tab$year$coef.glm),
              Elevation = exp(glmSRnat.overall$covar.tab$DEM_10$coef.glm),
              Slope = exp(glmSRnat.overall$covar.tab$SLOPE$coef.glm), 
              Northness = exp(glmSRnat.overall$covar.tab$Northern$coef.glm),
              SRali = exp(glmSRnat.overall$covar.tab$SRali$coef.glm), 
              FocalSpecies = exp(rowMeans(glmSRnat.overall$est))),
        outline = FALSE, ylab = "Effect size (exp(GLM coef))", cex.axis = 0.7)
abline(h = 1, col = "lightgrey")

# Frequency of significance for each factor:

barplot(cbind(  Year = table(!glmSRnat.overall$covar.tab$year$P.coef<0.05),
                Elevation = table(!glmSRnat.overall$covar.tab$DEM_10$P.coef<0.05),
              Slope =  table(!glmSRnat.overall$covar.tab$SLOPE$P.coef<0.05), 
               Northness =  table(!glmSRnat.overall$covar.tab$Northern$P.coef<0.05),
              SRali =  table(!glmSRnat.overall$covar.tab$SRali$P.coef<0.05),
              FocalSpecies = table(!rowSums(glmSRnat.overall$P<0.05, na.rm = T)>=1)),
        ylab = "Number of focal species", col = c("grey20","white"), cex.names = 0.7)






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







#______________________________________FIGURES for the paper ______________________________________###

#### Figure 1: barplot of frequencies:    #########

# extract frequency table:
effects <- glmSRnat.sum$class.summary
n <- length(unique(effects$group))

#draw barplot
par(mfrow = c(1,1), mar=c(5,3,3,1), las=1)
S <- effects[effects$group == "ALIEN:1",]
S <- rbind( c(NA,1,47,47, rep(0, dim(S)[2]-4)), S) ## Add the "rare" first class for illustration

barplot(S$nb.sp, ylim=c(0,50),col= "grey80",  border= NA, axes=F)
par(new=T)
barplot(S$freq.negative.above, ylim=c(0,50),col= "grey60",  border= NA, axes=F)
par(new=T)
b <- barplot(S$freq.thr, col="black" , ylim=c(0,50), border= NA, axes=F)
axis(2, tcl = -0.3, cex.axis = 0.8, mgp = c(1,0.5,0), las = 1)
text(y=-1.5, x = b+0.3, labels= abclasses[1:6],  cex=0.8, srt=45, adj=1, xpd = NA)
legend(x=0, y=60, bty="n", bg="white", 
       legend=c( "Total number of focal species",
                 "Number of negative effects",
                 "Number of critical abundances"),
       fill=c( "grey80","grey60","black"), border= c("grey90","grey60","black"), cex=0.7, xpd=NA,y.intersp =1)
mtext(text="number of species", side=2, outer=F, line=2, las=0, cex=1)
mtext(text="Abundance class", side=1, outer=F, line=3.5)




#### Figure 2: trends in alpha richness effect size ####

par(mfrow = c(2,4), mar=c(0,0,2,1), oma=c(7,7,1,1))
sel = impsp
M <- glmSRnat.overall
ylim=c(-100,50)
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
  col[(1:6)>= as.numeric( M$impact.spread[sp,threshold])] <- "black"  # above threshold
  col[(2:6)[hi>0]] <- "white"  # robust negative coef
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
  
  #threshold line
  # abline(v = M$impact.spread[sp,threshold]-1, lty="dashed")
  th <-M$impact.spread[sp,threshold]
  # abline(v = M$impact.spread[sp,threshold]-1, lty="dotted", lwd = 0.90, col = "grey60")
  arrows(x0 =th-1, y0 = 20,x1 =th-1, y1 = 35 , angle = 50, length = 0.05, col = "grey50", lwd = 2, code = 1, lend=0)
  
  # draw the points and lines
  par(new=T)
  plot(1:5,es, bg = col[2:6], pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F, bty = "n")
  
  
  # Add species name
  mtext(3, text=paste(letters[i],") ", species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.1, cex=0.7, line=0.2, las = 1)
  
  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  if (i %in% c(5:8)) text(y=-115, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)

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


par(mfrow = c(2,4), mar=c(0,0,2,1), oma=c(7,7,1,1))

ylim = c(-45,45)
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
  
  lines(1:5, lownull,col="grey70", lty="solid")
  lines(1:5, hinull,col="grey70", lty="solid")
  
  # if (!all(is.na(small))) points(1:6, small, pch=21, cex = 0.6)
  
  # annotations

  # X axis labels
  axis(1,at = 1:5, labels = F, tcl= 0.1,mgp=c(1,0.5,0),las=1, lwd = 0, lwd.ticks = 1)
  # if (i %in% c(8:11)) text(y=-53.2, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
  if (i %in% c(5:8)) text(y=-53.2, x = 1:5, labels= abclasses[2:6],  cex=1, srt=45, adj=1, xpd = NA)
 
  # Y axis
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, labels = F, lwd = 0, lwd.ticks = 1)
  if ( i %in% c(1,5)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1, lwd=0)
  
  # Y axis name
  if ( i %in% c(1,5)) {
      mtext(2, text=substitute(Delta*"("*gamma*"-richness)"[a], env = list(a = "rare")), adj=0.5, cex=0.8, line=2, las = 0)
  }
  
  
  
  # dotted horizontal
  abline(h=0,lty="dotted")
  
  #plot observed points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F, bty="n")
  box(bty="o", lwd = 1)
  
  # Add arrow indicating threshold for alpha richness
  th <- glmSRnat.overall$impact.spread[sp,threshold]
  # arrows(x0 =th-1, y0 = as.numeric(y[th] +20) ,x1 =th-1, y1 = as.numeric(y[th] +10), length = 0.07, col = "grey60", lwd = 2, code = 2)
  # arrows(x0 =th-1, y0 = -35,x1 =th-1, y1 = -42, length = 0.07, col = "grey60", lwd = 1, code = 1)
  # abline(v = M$impact.spread[sp,threshold]-1, lty="dotted", lwd = 0.90, col = "grey60")
  # arrows(x0 =th-1, y0 = 25,x1 =th-1, y1 = 30 , angle = 50, length = 0.05, col = "grey60", lwd = 3, code = 1, lend=0)
  arrows(x0 =th-1, y0 = 30,x1 =th-1, y1 = 40 , angle = 50, length = 0.05, col = "grey50", lwd = 2, code = 1, lend=0)
  
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


#### Figure 4: NOT RUN Maps of presence and impact  WITH LAND AND OCEAN ####

# library (prettymapr)
# db = databp[databp$PlotName %in% unimprovedgrasslands, ]
# targets <- impsp
# region <- study_area
# par(
#   mfrow = c(3,3),
#   mar = c(0, 0, 1, 0),
#   oma = c(2, 0, 0, 0)
# )
# 
# 
# 
# for (i in 1:9) {
#   if (i ==9) {
#     plot(region, border = NA)
#     legend(
#       'left',
#       legend = c("< critical abundance", "> critical abundance", "dominant"),
#       fill = c("tan", "sienna2", "sienna4"),
#       xpd = NA,
#       bty = "n",
#       border = c("tan", "sienna2", "sienna4")
#     ) 
#   
#     addnortharrow(pos = "topright", padin = c(0.40, 0.4), scale = 0.5,
#                              lwd = 1, border = "black", cols = c("white", "black"),
#                              text.col = "black")
# 
#     map.scale(x = extent(region)[1]+0.1, y = extent(region)[3] + 0.05, ratio = F, relwidth = 0.40) 
#   }
#   
#   if (i <9) {
#     sp=targets[i]
#     plot.alien = as.character(db[which(db$abun %in% c(1, 2, 3, 4, 5, 6) &
#                                          db$vegtype == "G" & db$SpeciesCode == sp), "PlotName"])
#     th = glmSRnat.overall$impact.spread[sp, "th.CI"]
#     es <- as.numeric(glmSRnat.overall$est[sp, ])
#     n <- as.numeric(glmSRnat.overall$n.obs[sp, ])[2:6]
#     es[which(n < 5)] <- NA
#     es = which(!is.na(es)) + 1
#     m <- max(es, na.rm = T)
#     max.impact = as.character(db[which(db$abun >= 6 &
#                                          db$vegtype == "G" & db$SpeciesCode == sp), "PlotName"])
#     th.impact = as.character(db[which(db$abun >= th &
#                                         db$vegtype == "G" & db$SpeciesCode == sp), "PlotName"])
#      
#     
#     plot(region, col = "grey50", border = NA)
#     polygon(extent(region)[c(1,2,2,1)], extent(region)[c(3,3,4,4)], col ="lightblue", border = NA)
#     plot(region, add= T, col = "grey70", border = NA)
#     plot(study_area, add= T, col = "grey90", border = NA)
#     plot(region, add= T, border = "grey70")
#     
#     plot(
#       envplot[plot.alien,],
#       pch = 22,
#       cex = 0.5,
#       col = "tan",
#       bg = "tan",
#       add = T
#     )
#     plot(
#       envplot[th.impact,],
#       pch = 22,
#       cex = 0.5,
#       col = "sienna2",
#       bg = "sienna2",
#       add = T
#     )
#     plot(
#       envplot[max.impact,],
#       pch = 22,
#       cex = 0.5,
#       col = "sienna4",
#       bg = "sienna4",
#       add = T
#     )
#     
#     mtext(
#       3,
#       text = paste(letters[i],") ",sub("_", " ", species[sp, "tip"]), sep = ""),
#       adj = 0.1,
#       line = -0.5,
#       font = 3,
#       cex = 0.6
#     )
#   }
#     
# }
# 
# 


#### Figure 5: correlation of deltagamma with delta alpha and with spread  ######

#loss in alpha vs. loss in gamma (lg transformed pearson corr)   
da = -da
dg = -dg
par(mfrow = c(1,2), mar=c(3,1,2,1), oma=c(1,2,0,1), las = 1, xpd = TRUE)
plot(da, dg ,pch = 20, col=cgam , ann=F,type ="n", axes =F, log = '', xlim = c(0,4))
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
text(x = da,y = dg,rownames(table.div.part) , col="grey", cex =0.7, pos = 4)
points(x = da,y = dg, pch = 20)
box(bty="l")
f <- cor.test(da ,dg, method = "spearman")
mtext(3, text = substitute(italic(rho == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=F))),adj = 1, cex=0.8, font = 3, line= 0.5)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1, line= 0.5)
mtext(1,text = expression(Delta*alpha*"-richness"), line=2 )
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=2 , las=0)

### spread above critical abundance - log transforming the spread of potential impact
print(f <- cor.test(spread ,dg, method = "spearman"))

plot (spread,dg , pch = 20, col=cgam , ann=F, log = "", xlim =c(5, 500), type ="p", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=F, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
dg.lab <- dg
text(x = spread + 12, y = dg.lab - c(0,2,0,0,0,0,0,0) ,label = rownames(table.div.part) , col="grey", cex =0.7,  adj = c(0,0))
points(x = spread, y = dg, pch = 20)
box(bty="l")
mtext(3, text = substitute(italic(rho == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8, font = 3, line= 0.5)
mtext(1,text = "Number of plots > critical abundance", line=2 )
mtext(3, text = 'b)',adj = 0, cex=0.8, font = 1, line= 0.5)
abline(lm(dg ~ spread), xpd = FALSE, lty = "dotted")
  
  
  
#### Additional regressions ####
### spread above critical abundance - log transforming the spread of potential impact
print(f <- cor.test(dbeta ,dg, method = "spearman"))

plot (dbeta,dg , pch = 20, col=cgam , ann=F,  type ="p", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=F, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
dg.lab <- dg
text(x = dbeta, y = dg.lab ,label = rownames(table.div.part) , col="grey", cex =0.6, pos = 4)
points(x = dbeta, y = dg, pch = 20)
box(bty="l")
mtext(3, text = substitute(italic(rho == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8, font = 3, line= 0.5)
mtext(1,text = expression(Delta*beta*"turnover"["c"]), line=2 )
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=2 , las=0)
mtext(3, text = 'c)',adj = 0, cex=0.8, font = 1, line= 0.5)



# Above.Acrit*beta is the best:
cor.test( formula =  ~ delta.gamma.c + above.Acrit, data = mat, method = "spearman")
cor.test( formula =  ~ delta.gamma.c + Presence, data = mat, method = "spearman")
cor.test( formula =  ~ delta.gamma.c + dominance, data = mat, method = "spearman")


regsubsets(formula =  delta.gamma.c ~ above.Acrit * delta.beta.z, data = mat)
# Changes in -richness above critical abundances were best explained (exhaustive variable selection using AIC) by a model including the interaction between the decrease in -dissimilarity and the number of plots above critical abundance: loss in -richness was highest when a large number of plots above critical abundance were combined with a larger loss in -dissimilarity (Δc ~ Nb.Plots(>Acrit) *  Δc ; R2  = 0.83, df = 4, P  < 0.05).


summary(f0 <- lm( formula =  delta.gamma.c ~  above.Acrit , data = mat))
summary(f0b <- lm( formula =  delta.gamma.c ~  delta.beta.z , data = mat))
summary(f1 <- lm( formula =  delta.gamma.c ~  above.Acrit + delta.beta.z , data = mat))
summary(f <- lm( formula =  delta.gamma.c ~ above.Acrit * delta.beta.z, data = mat))

AIC(f0, f0b, f1, f)
BIC(f0, f0b, f1, f)
anova(f)
as.matrix(anova(f))[,2] /sum(as.matrix(anova(f)[,2]))  # % sum of squares explained

plot(-dg, predict(f), ylim= c(-60,0), xlim= c(-60,0))
abline(0,1)

tmp <- mat
tmp$delta.gamma.c <- - tmp$delta.gamma.c
tmp$delta.beta.z <- - tmp$delta.beta.z
summary(f <- lm( formula =  delta.gamma.c ~ above.Acrit * delta.beta.z, data = tmp))
as.matrix(anova(f))[,2] /sum(as.matrix(anova(f)[,2]))  # % sum of squares explained

# Prediction based on correlation between predictors:

summary(fb <- lm(delta.beta.z ~ above.Acrit, mat ))
tmp <- data.frame(above.Acrit = 1:500, delta.beta.z = coef(fb)[2] * 1:500 + coef(fb)[1])
summary(f <- lm( formula =  delta.gamma.c ~  above.Acrit * delta.beta.z, data = mat))
pr <- predict(f, newdata= tmp, se.fit = T)

par(mfrow=c(1,1), xpd = FALSE)
plot (spread,-dg , pch = 20, col=cgam , ann=F, log = "", xlim =c(5, 500), type ="p", axes =F)
axis(1, mgp=c(0,0.3,0), tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)
axis(2, mgp=c(0,0.3,0),labels=T, tcl=0.2 ,cex.axis = 0.9, lwd = 0, lwd.ticks = 1)

points(x = spread, y = dg, pch = 20)
box(bty="l")
mtext(3, text = substitute(italic(rho == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=T))),adj = 1, cex=0.8, font = 3, line= 0.5)
lines(tmp$above.Acrit, pr$fit, col = "grey")
lines(tmp$above.Acrit, pr$fit + pr$se.fit, col = "grey", lty = "dotted")
lines(tmp$above.Acrit, pr$fit - pr$se.fit,  col = "grey", lty = "dotted")

dg.lab <- dg
text(x = spread + 12, y = dg.lab - c(0,2,0,0,0,0,0,0) ,label = rownames(table.div.part) , col="darkgrey", cex =0.6,  adj = c(0,0))

mtext(1,text = "Number of plots > critical abundance", line=2 )
mtext(3, text = 'b)',adj = 0, cex=0.8, font = 1, line= 0.5)
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=1.5 , las=0)

## Plot da and da.max

plot(da, da.max)
text(da, da.max, labels = names(da.max))
abline(0,1)
cor.test(da, da.max, method = "spearman")
plot(dg, da.max - da)
plot(critical.abun, da.max)
