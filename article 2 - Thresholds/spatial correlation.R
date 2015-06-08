#### Geographical clustering/spatial autocorrelation of impacted plots

# geographical distance matrix
spDst <-pointDistance(coordinates(envplot), lonlat = F, allpairs=FALSE)
rownames(spDst) <- colnames(spDst) <-rownames(coordinates(envplot))


# autocorrelation among impacted and non impacted sites


## home made null model of spatial mean nearest distance
spatial.mnnd <- function(sp.dist = spDst, X = var==2, nreps = 999) {
  stopifnot(dim(sp.dist)[1] == length(X))
diag(sp.dist) <- NA
mnnd.obs <- mean(apply(sp.dist[X, X], 1, min, na.rm=T))

mnnd.rand = rep(NA, nreps)
for (i in 1:nreps){
  rv <- sample(X)
  mnnd.rand[i] <- mean(apply(sp.dist[rv, rv], 1, min, na.rm=T))
}

mnnd.rand <-  c(mnnd.obs,  mnnd.rand)
hist(mnnd.rand)
abline(v = mnnd.obs, col ="red")

p.mnnd <- (sum(mnnd.rand <mnnd.obs) + sum(mnnd.rand ==mnnd.obs)/2)/(nreps+1)


return(cbind(obs = mnnd.obs, mean.null = mean(mnnd.rand), sd.null = sd(mnnd.rand),P = p.mnnd))
}



# apply to target species threshold

mnnd.th <- matrix(NA, nrow= 11, ncol = 8)
colnames(mnnd.th) <- paste(c("mnnd.obs","null.mean", "null.sd", "P<obs"), 
               c(rep("below", 4), rep("above", 4)), sep = "_")
rownames(mnnd.th) <- impsp

for (i in 1:length(impsp)) {

sp=impsp[i]


community <- comm
 community <- community[!is.na(community[,sp]) & community[,sp]>0,]

dst <- spDst[rownames(community),rownames(community)]

var <- community[,sp]
var[var < glmSRnat.overall$impact.spread [sp,"th.CI"] & var!=0] <- 1
var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 2


mnnd.th [i,1:4] <- spatial.mnnd(sp.dist = dst, X = var==1, nreps = 999)
mnnd.th [i,5:8] <- spatial.mnnd(sp.dist = dst, X = var==2, nreps = 999)
print(paste(i, sp))
}









#geary test

library(ade4)
randGeary <- gearymoran(spDst, X = community[, impsp], nrepet =999, alter=c("two-sided"))


# varigram
breaks = seq(0, max(spDst), l = 100)
v1 <- variog(coords =coord(envplot)[var==2],data = envplot$SRnat,breaks = breaks)


