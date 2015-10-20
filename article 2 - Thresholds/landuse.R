## Land use effects on trends

### plot proportion of low vs high productivity grassland in each abundance class for 11 species
abclasses= c("Rare" ,"Occasional",  "Frequent", "Common", "Abundant","Dominant")
threshold ="th.CI"
sel <- impsp
M <- glmSRnat.overall

### proportion of low productivity grasslands in each abundance class
par(mfrow = c(3,4), mar=c(2,2,2,2), oma=c(3,3,2,3))
for (i in 1:length(sel))  {
  sp <- sel[i]

ab.vector <- data.frame(abun = comm[realgrasslands,sp])
rownames(ab.vector) <- rownames(comm[realgrasslands,])
ab.vector$landuse <- envplot@data[rownames(ab.vector),]$landuse
ab.vector$SRnat <- envplot@data[rownames(ab.vector),]$SRnat
ab.vector$ALIEN.dom  <- envplot@data[rownames(ab.vector),]$ALIEN.dom
ab.vector$grasslands <- "other"
ab.vector$grasslands[grep("Grassland - With woody biomass",ab.vector$landuse)] <- "woody"
ab.vector$grasslands[grep("High",ab.vector$landuse)] <- "high"
ab.vector$grasslands[grep("Low",ab.vector$landuse)] <- "low"

tab <- t(as.data.frame.matrix(table(ab.vector$landuse, ab.vector$abun)))

tmp <- barplot(t(tab/rowSums(tab)), col=LUCAS.col, main=sp,    cex.axis=0.8, cex.names=0.8, las=1)
abline(v =tmp[glmSRnat.overall$impact.spread[sp,threshold]+1], lty="dashed")
}
mtext(2, outer=T, text = "proportion of improved pastures", line =1)


### proportion of individual species incidence in each class
par(mfrow = c(3,4), mar=c(2,2,2,2), oma=c(3,3,2,3))
for (i in 1:length(sel))  {
  sp <- sel[i]

  ab.vector <- data.frame(abun = comm[realgrasslands,sp])
  rownames(ab.vector) <- rownames(comm[realgrasslands,])
  ab.vector$landuse <- envplot@data[rownames(ab.vector),]$landuse
  ab.vector$SRnat <- envplot@data[rownames(ab.vector),]$SRnat
  ab.vector$ALIEN.dom  <- envplot@data[rownames(ab.vector),]$ALIEN.dom
  ab.vector$grasslands <- "other"
  ab.vector$grasslands[grep("High",ab.vector$landuse)] <- "high"
  ab.vector$grasslands[grep("Low",ab.vector$landuse)] <- "low"


  ab.vector$first.rank  <- envplot@data[rownames(ab.vector),]$first.rank
  ab.vector$lolper.pres <- comm[realgrasslands,"LOLPER"]>0
  ab.vector$trirep.pres <- comm[realgrasslands,"TRIREP"]>0
  ab.vector$achmil.pres<- comm[realgrasslands,"ACHMIL"]>0

  tab <- table(ab.vector$grasslands, ab.vector$abun)
  plot(colnames(tab), tab[2,] / colSums(tab),pch=23, bg = "green", col = "green", main = sp, ylim = c(0,1))

  tab3 <- table(ab.vector$lolper.pres, ab.vector$abun)
  points(colnames(tab), tab3[2,] / colSums(tab),pch=20, cex=1.5, col = "blue", main = sp, ylim = c(0,1))
  tab4 <- table(ab.vector$trirep.pres, ab.vector$abun)
  points(colnames(tab), tab4[2,] / colSums(tab),pch=20, cex=1.5, col = "brown", main = sp, ylim = c(0,1))
  tab4 <- table(ab.vector$achmil.pres, ab.vector$abun)
  points(colnames(tab), tab4[2,] / colSums(tab),pch=20, cex=1.5, col = "pink", main = sp, ylim = c(0,1))

  abline(v =glmSRnat.overall$impact.spread[sp,threshold], lty="dashed")
}


### correlation between abundance of the 11 impsp

cor.impsp <- cor(comm[, impsp], use = "complete.obs", method= "spearman")
heatmap(cor.impsp)

heatmap(ceiling(cor.impsp>0.2))

heatmap(as.matrix(dist(t(comm[, impsp]))))
heatmap(ceiling(t(comm[, impsp] >0))
