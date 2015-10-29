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



### test proportion of improved grassland at rare vs. threshold abundance
par(mfrow = c(3,4), mar=c(2,2,2,2), oma=c(3,3,2,3))

pasture.effect<- matrix(NA, nrow= 12, ncol = 10)

for (i in 1:length(sel))  {
  sp <- sel[i]

  ab.vector <- data.frame(abun = comm[realgrasslands,sp])
  rownames(ab.vector) <- rownames(comm[realgrasslands,])
  ab.vector$landuse <- envplot@data[rownames(ab.vector),]$landuse
  ab.vector$SRnat <- envplot@data[rownames(ab.vector),]$SRnat
  ab.vector$ALIEN.dom  <- envplot@data[rownames(ab.vector),]$ALIEN.dom
  ab.vector$grasslands <- "other"
  ab.vector$grasslands[grep("High",ab.vector$landuse)] <- "high"


  tab <- t(as.data.frame.matrix(table(ab.vector$grasslands, ab.vector$abun)))
  th<- glmSRnat.overall$impact.spread[sp,threshold]

  prop.tab <- as.matrix(tab/rowSums(tab))
tbl <-  table(ab.vector$grasslands, ab.vector$abun)[,c(2,th+1)]
  (cor.pasture <- cor.test(prop.tab[2:dim(tab)[1],1],2:dim(tab)[1], method = "spearman", exact = FALSE))
  (t.pasture <- t.test(tbl))
  (Ftest.pasture <-  fisher.test(tbl,simulate.p.value = TRUE))



  pasture.effect[i,] <- c(th, as.numeric(prop.tab[c(2,th+1),1]),
    S.rho = cor.pasture$est, S.P = cor.pasture$p.val,
                                     t =t.pasture$statistic, t.df = t.pasture$parameter, t.P = t.pasture$p.value,
                                     F.odds =Ftest.pasture$estimate, F.P = Ftest.pasture$p.value)


  tmp <- barplot(t(tab/rowSums(tab)), col=c("grey", "white"), main=sp,    cex.axis=0.8, cex.names=0.8, las=1)
  abline(v =tmp[glmSRnat.overall$impact.spread[sp,threshold]+1], lty="dashed")

  }
rownames(pasture.effect) <- impsp

colnames(pasture.effect) <- c("th","p.rare","p.th","S.rho", "S.P","t", "t.df", "t.P", "F.odds", "F.P")

