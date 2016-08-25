### gamma rarefaction curves
# OUTDATED see gamma estimators for updated version

library(vegan)
### my rarefy function : rarefies the number of plots sampled
rarefy.gamma <- function (x, sample , nreps, replace =T) {
  g <- rep(NA,nreps)
  nr <- nrow(x)
for ( i in 1:nreps) {
 spl <- sample(1:nr, sample, replace = T)
 g[i] <- sum(colSums(x[spl,])>0, na.rm =T)
}
  return(g)
}

gamma.rarecurve <- function (com, step=20, nreps = 19, plot.curve =F, ...){

n.row <- nrow(com)
sample.vector <- seq(2, n.row, by = step)
if (max(sample.vector) < n.row) sample.vector <- c(sample.vector, n.row)

gamma.curve <- matrix(NA, nrow = nreps, ncol = length(sample.vector))

for (s in 1:length(sample.vector)) {
  gamma.curve[,s] <- rarefy.gamma(x = com, sample = sample.vector[s], nreps = nreps,replace =replace)
  print(paste("gamma curve size:",s, "of", length(sample.vector) ))
}

curve.summary <- data.frame( samples = sample.vector,
                             mean.gammas = colMeans(gamma.curve),
                             sd.gammas = apply(gamma.curve, 2, sd, na.rm=T))

if (plot.curve) plot(mean.gammas ~ samples, curve.summary, type ="l", ...)
 return(curve.summary)
}

gamma.rarecurve (community)

### Rarefaction curves above / below
group=natives
par(mfrow = c(3,4), mar=c(2,2,2,2), oma=c(4,4,2,2))
for(sp in impsp) {
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,colSums(community)>0]
  community<- community[,colnames(community) %in% group]
  community<- ceiling(community>0)

  var <- comm[rownames(community),sp]
  abun <- table(var)
  abun <- as.numeric(names(abun[abun>=5]))


  th <- glmSRnat.overall$impact.spread[sp,]$th.CI
  th.var <- rep(0, length(var))
  th.var[var>=th] <-1

  curve.below <- gamma.rarecurve(com = community[th.var == 0,], step = 5, nreps =29, plot.curve=F, replace = F)
  curve.above <- gamma.rarecurve(com = community[th.var == 1,], step = 5, nreps =29, plot.curve=F, replace = F)


  # rarefaction curves for total gamma richenss in each class
  plot(c(1,max(curve.above$samples, curve.below$samples)),
       c(1,max(curve.above$mean.gammas, curve.below$mean.gammas)), type="n", ann=F)
  lines(spline(curve.below))
  points(curve.below)
  lines(spline(curve.above), col="red")
  points(curve.above)
   mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)

}

### Rarefaction curve for one species abundance gradient

par(mfrow = c(3,4), mar=c(2,2,2,2), oma=c(4,4,2,2))
for(sp in impsp) {
community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
community <- community[,colSums(community)>0]
community<- community[,colnames(community) %in% group]
community<- ceiling(community>0)

var <- comm[rownames(community),sp]
abun <- table(var)
abun <- as.numeric(names(abun[abun>=5]))

# create community matrix of sp. frequency in each abundance class
freq <- colSums(community) / dim(community)[1]
freq.ab <- sapply(1:6,function(k) {
 f <- (if (sum(var==k)>0) sapply(colnames(community), function(s) sum(community[which(var==k), s], na.rm=T))
else rep(NA, length(colnames(community))) )
  })

# number of plots in each class
n.plots <- sapply(1:6,function(k) {
 sum(var==k)
})

# variation in species frequencies
# plot( c(1,6),c(0,1), type= "n")
# sapply(rownames(freq.ab), function(spi){
#  lines( 1:6, freq.ab[spi,]/n.plots)
# })

# rarefaction curves for total gamma richenss in each class
com.freq <- t(freq.ab)
rownames(com.freq) <- abclasses
com.freq <-com.freq[which(rowSums(com.freq)>5),]

rarecurve(com.freq, step = 2)
mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
      font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)

}]


par(mfrow = c(3,4), mar=c(2,2,2,2), oma=c(4,4,2,2))
for(sp in impsp) {
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,colSums(community)>0]
  community<- community[,colnames(community) %in% group]
  community<- ceiling(community>0)

  var <- comm[rownames(community),sp]
  abun <- table(var)
  abun <- as.numeric(names(abun[abun>=5]))


  th <- glmSRnat.overall$impact.spread[sp,]$th.CI
  th.var <- rep(0, length(var))
  th.var[var>=th] <-1
  # create community matrix of sp. frequency in each abundance class
  freq <- colSums(community) / dim(community)[1]

  freq.ab <- sapply(c(0,1),function(k) {
    f <- (if (sum(th.var)>0) sapply(colnames(community), function(s) sum(community[which(th.var==k), s], na.rm=T))
          else rep(NA, length(colnames(community))) )
  })

  # number of plots in each class
  n.plots <- table(th.var)

  # variation in species frequencies
  # plot( c(1,6),c(0,1), type= "n")
  # sapply(rownames(freq.ab), function(spi){
  #  lines( 1:6, freq.ab[spi,]/n.plots)
  # })

  # rarefaction curves for total gamma richenss in each class
  com.freq <- t(freq.ab)
  rownames(com.freq) <- c("below", "above")
  com.freq <-com.freq[which(rowSums(com.freq)>5),]

  rarecurve(com.freq, step = 2, col=c("black", "red"))
  mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)

}
