### Species covariations

### proportion of individual species incidence in each class
par(mfrow = c(3,4), mar=c(2,2,2,2), oma=c(3,3,2,3))
for (i in 1:length(sel))  {
  sp <- sel[i]

  ab.vector <- data.frame(abun = comm[unimprovedgrasslands,sp])
  rownames(ab.vector) <- rownames(comm[unimprovedgrasslands,])
  ab.vector$landuse <- envplot@data[rownames(ab.vector),]$landuse
  ab.vector$SRnat <- envplot@data[rownames(ab.vector),]$SRnat
  ab.vector$ALIEN.dom  <- envplot@data[rownames(ab.vector),]$ALIEN.dom
  ab.vector$grasslands <- "other"
  ab.vector$grasslands[grep("High",ab.vector$landuse)] <- "high"
  ab.vector$grasslands[grep("Low",ab.vector$landuse)] <- "low"


  ab.vector$first.rank  <- envplot@data[rownames(ab.vector),]$first.rank
  ab.vector$lolper.pres <- comm[unimprovedgrasslands,"LOLPER"]>0
  ab.vector$trirep.pres <- comm[unimprovedgrasslands,"TRIREP"]>0
  ab.vector$achmil.pres<- comm[unimprovedgrasslands,"ACHMIL"]>0

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


#### Plot abundance correlation between all 18 impsp
n <- unique(c(impsp,impsp.ali))
tmp <- comm[unimprovedgrasslands,n]
tmp <- tmp[-which(rowSums(tmp)==0),]
tmp[tmp ==0] <- NA   # remove absence information

cortest.allimp<- list( r = matrix(NA, nrow=18, ncol = 18, dimnames = list(n, n)),
                          df = matrix(NA, nrow=18, ncol = 18, dimnames = list(n, n)),
                          P = matrix(NA, nrow=18, ncol = 18, dimnames = list(n, n)))
par(mfrow=c(18,18), mar=c(1,1,1,1), oma=c(0,6,2,0))
for (i in 1:length(n)){
  spi <- n[i]
  for (j in 1:length(n)) {
    spj <- n[j]
    y= tmp[,spi]
    x= tmp[,spj]
    common = !(is.na(x)) & !(is.na(y))
    x<- x[common]
    y<- y[common]
    if (sum(common)<3)  cortest.allimp$r[i,j] <- 0
    if (sum(common)>=3) {
      ct <- cor.test(x,y,use = "pairwise.complete.obs", method = "pearson", exact=FALSE)
      cortest.allimp$r[i,j] = ct$est
      cortest.allimp$P[i,j] = ct$p.value
      cortest.allimp$df[i,j] = ct$parameter
    }

#     plot(jitter(x), jitter(y), ann=F, axes=F, col="#ababab50", bg="#ababab50",pch=21)
#     if(ct$p.value <=0.05 & !is.na(ct$p.value)) abline(lm(y~x))
#     mtext(3, text = paste(round(ct$est,2), p2star(ct$p.value)), cex= 0.7, adj=1)
#
#     if (j==1) mtext(2, text=spi, line = 1, las=1, cex= 0.7)
#     if (i==1) mtext(3, text=spj, line = 1, las=0, cex=0.7)

  }
}

out <- sapply(n, FUN = function(sp) {
  paste(round(cortest.allimp$r[sp,],2),sapply(cortest.allimp$P[,sp], p2star), sep ="")
})
rownames(out) <- n

write.csv(out,file = "corr.allimp.csv")




#### Plot abundance correlation between 12 species from impsp
tmp <- comm[unimprovedgrasslands, impsp]
tmp <- tmp[-which(rowSums(tmp)==0),]
tmp[tmp ==0] <- NA   # remove absence information


cortest.impsp<- list( r = matrix(NA, nrow=12, ncol = 12, dimnames = list(impsp, impsp)),
                      df = matrix(NA, nrow=12, ncol = 12, dimnames = list(impsp, impsp)),
                      P = matrix(NA, nrow=12, ncol = 12, dimnames = list(impsp, impsp)))
par(mfrow=c(12,12), mar=c(1,1,1,1), oma=c(0,6,2,0))
for (i in 1:length(impsp)){
  spi <- impsp[i]
  for (j in 1:length(impsp)) {
    spj <- impsp[j]
    y=tmp[,spi]
    x= tmp[,spj]
    common = !(is.na(x)) & !(is.na(y))
    x<- x[common]
    y<- y[common]

    ct <- cor.test(x,y,use = "pairwise.complete.obs", method = "pearson", exact=FALSE)
    cortest.impsp$r[i,j] = ct$est
    cortest.impsp$P[i,j] = ct$p.value
    cortest.impsp$df[i,j] = ct$parameter
    plot(jitter(x), jitter(y), ann=F, axes=F, col="#ababab50", bg="#ababab50",pch=21)
    if(ct$p.value <=0.05) abline(lm(y~x))
    mtext(3, text = paste(round(ct$est,2), p2star(ct$p.value)), cex= 0.7, adj=1)

    if (j==1) mtext(2, text=spi, line = 1, las=1, cex= 0.7)
    if (i==1) mtext(3, text=spj, line = 1, las=0, cex=0.7)
  }
}

out <- sapply(impsp, FUN = function(sp) {
  paste(round(cortest.impsp$r[sp,],2),sapply(cortest.impsp$P[,sp], p2star), sep ="")
})
rownames(out) <- impsp

write.csv(out,file = "corr.impsp.csv")

#### Plot abundance correlation between 10 species from impsp.ali
n <- impsp.ali
tmp <- comm[unimprovedgrasslands,n]
tmp <- tmp[-which(rowSums(tmp)==0),]
tmp[tmp ==0] <- NA   # remove absence information

cortest.impsp.ali<- list( r = matrix(NA, nrow=10, ncol = 10, dimnames = list(n, n)),
                      df = matrix(NA, nrow=10, ncol = 10, dimnames = list(n, n)),
                      P = matrix(NA, nrow=10, ncol = 10, dimnames = list(n, n)))
par(mfrow=c(10,10), mar=c(1,1,1,1), oma=c(0,6,2,0))
for (i in 1:length(n)){
  spi <- n[i]
  for (j in 1:length(n)) {
    spj <- n[j]
    y= tmp[,spi]
    x= tmp[,spj]
    common = !(is.na(x)) & !(is.na(y))
    x<- x[common]
    y<- y[common]
    if (sum(common)<3)  cortest.allimp$r[i,j] <- 0
    if (sum(common)>=3) {
    ct <- cor.test(x,y,use = "pairwise.complete.obs", method = "pearson", exact=FALSE)
    cortest.impsp.ali$r[i,j] = ct$est
    cortest.impsp.ali$P[i,j] = ct$p.value
    cortest.impsp.ali$df[i,j] = ct$parameter
    }

    plot(jitter(x), jitter(y), ann=F, axes=F, col="#ababab50", bg="#ababab50",pch=21)
    if(ct$p.value <=0.05 & !is.na(ct$p.value)) abline(lm(y~x))
    mtext(3, text = paste(round(ct$est,2), p2star(ct$p.value)), cex= 0.7, adj=1)

    if (j==1) mtext(2, text=spi, line = 1, las=1, cex= 0.7)
    if (i==1) mtext(3, text=spj, line = 1, las=0, cex=0.7)

  }
}

out <- sapply(impsp.ali, FUN = function(sp) {
  paste(round(cortest.impsp.ali$r[sp,],2),sapply(cortest.impsp.ali$P[,sp], p2star), sep ="")
})
rownames(out) <- impsp.ali

write.csv(out,file = "corr.impsp.ali.csv")



d<- as.dist(1-cortest.impsp.ali$r)
d[is.na(d)] <- 0
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=4) # cut tree into 5 clusters
# draw dendogram with red borders around the 4clusters
rect.hclust(fit, k=4, border="red")




## cluster analysis
# Ward Hierarchical Clustering
tmp <- comm[unimprovedgrasslands, impsp]
tmp <- tmp[-which(rowSums(tmp)==0),]
tmp[tmp ==0] <- NA   # remove absence information

d <- vegdist(t(tmp), method = "bray",na.rm = TRUE) # distance matrix

fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=4) # cut tree into 5 clusters
# draw dendogram with red borders around the 4clusters
rect.hclust(fit, k=4, border="red")


# k means
tmp <- comm[unimprovedgrasslands, impsp]
tmp <- tmp[-which(rowSums(tmp)==0),]

# tmp[tmp ==0] <- NA   # remove absence information
fit <- kmeans(t(tmp), centers=4, iter.max = 100, nstart=25)
fit


# Determine number of clusters with sum of squares
mydata=t(tmp)
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:6) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:6, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

fit <- kmeans(mydata, 5) # 5 cluster solution
# get cluster means
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster)


# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)

tmp <- comm[unimprovedgrasslands, impsp]
tmp <- tmp[-which(rowSums(tmp)==0),]
tmp[tmp ==0] <- NA   # remove absence information
mydata = tmp  # this funtion clusters columns, not rows)
hclust.cor <- pvclust(mydata, method.hclust="ward",
               method.dist="correlation" )
plot(hclust.cor ) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(hclust.cor , alpha=.95)

## impsp ALIEN richness
tmp <- comm[unimprovedgrasslands, impsp.ali]
tmp <- tmp[-which(rowSums(tmp)==0),]

tmp[tmp ==0] <- NA   # remove absence information
mydata = tmp  # this funtion clusters columns, not rows)
hclust.cor.SRali <- pvclust(mydata, method.hclust="ward",
                      method.dist="correlation",  use.cor="pairwise.complete.obs")
plot(hclust.cor.SRali ) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(hclust.cor.SRali , alpha=.89)



diag(cortest.allimp$r)<-NA
heatmap(cortest.allimp$r, Rowv = as.dendrogram(hclust.cor$hclust),
        Colv =as.dendrogram(hclust.cor$hclust),
        col = colorRampPalette(c("slateblue","white", "firebrick"))(20))  ### uses the Ward.D

# coding positive and negative cassociations :
cor.impsp.simple <- ceiling(cortest.impsp$rho >0.15 &  cortest.impsp$P <=0.05) - ceiling(cortest.impsp$rho <0.15 &  cortest.impsp$P <=0.05)

heatmap(cor.impsp.simple, Rowv = as.dendrogram(hclust.cor$hclust),
        Colv =as.dendrogram(hclust.cor$hclust),
        col =colorRampPalette(c("slateblue","white", "firebrick"))(3))  ### uses the Ward.D
image(cor.impsp.simple,col =colorRampPalette(c("slateblue","white", "firebrick"))(3))


heatmap(cor.impsp.simple,col =colorRampPalette(c("slateblue","white", "firebrick"))(3))

# NMDS ordination of 12 species abundances  => inconclusive
# tmp[tmp==0] <- NA
# Ordinates mainly based on how frwquent they are (widespread species cluster together)= not good
## presence absence only
tmp <- comm[unimprovedgrasslands, impsp]
tmp <- tmp[-which(rowSums(tmp)==0),]
tmp[tmp>0] <- 1
nmds.impsp.sites.bin <- metaMDS(tmp, distance="jacc", k =3, trymax = 15, previous.best =nmds.impsp.sites.bin  )
nmds.k.species <- kmeans(nmds.impsp.sites.bin$species[, 1:3], centers= 4)
nmds.k.sites <- kmeans(nmds.impsp.sites.bin$points[, 1:3], centers= 4)

# with abundnances
tmp <- comm[unimprovedgrasslands, impsp]
tmp <- tmp[-which(rowSums(tmp)==0),]
nmds.impsp.sites.ab <- metaMDS(tmp, distance="bray", k =3, trymax = 15, previous.best =nmds.impsp.sites.ab)
nmds.k.ab.species <- kmeans(nmds.impsp.sites.ab$species[, 1:3], centers= 4,nstart = 4,iter.max=100)
nmds.k.ab.sites <- kmeans(nmds.impsp.sites.ab$points[, 1:3], centers= 4,nstart = 4, iter.max=100)

x11()
ordifit <- nmds.impsp.sites.bin
nmds.k.species <- kmeans(ordifit$species[, 1:3], centers= 4,nstart = 4,iter.max=100)
nmds.k.sites <- kmeans(ordifit$points[, 1:3], centers= 4,nstart = 4, iter.max=100)

lim=c(min(ordifit$species,ordifit$species)+0.2, max(nmds.impsp$points,nmds.impsp$species)+0.2)
plot(ordifit$points, col="darkgrey")
x <- ordifit$species[,1]
y <- ordifit$species[,2]
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col=nmds.k.species$cluster)
text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7, col =nmds.k.species$cluster )

lim=c(min(ordifit$species,ordifit$species)+0.2, max(nmds.impsp$points,nmds.impsp$species)+0.2)
plot(ordifit$points, col=nmds.k.sites$cluster)
x <- ordifit$species[,1]
y <- ordifit$species[,2]
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7,col="black" )
legend("topright", fill= unique(nmds.k.sites$cluster), legend=unique(nmds.k.sites$cluster))

boxplot( envplot@data[names(nmds.k.sites$cluster), "SRnat"] ~ nmds.k.sites$cluster, col=1:4)

### NMDS analysis of the entire invaded communities

# select realgrassland plots where at least one of the important target species is present:
community <- comm[which((rownames(comm) %in% unimprovedgrasslands) & (comm[,sp]>0) ),impsp]
community <- comm[which((rownames(comm) %in% unimprovedgrasslands) & (rowSums(comm[,impsp])>0) ),impsp]

nmds.results <- metaMDS(community, distance = "bray", k = 2, trymax = 10)

lim=c(min(nmds.results$points,nmds.results$species)+0.2, max(nmds.results$points,nmds.results$species)+0.2)
plot(nmds.results$points, col="darkgrey",ylim=lim, xlim=lim )
x <- nmds.results$species[,1]
y <- nmds.results$species[,2]
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7 )
i=1
arrows(rep(0), rep(0), x[i], y[i], length =0.1, col="red", cex=2)
