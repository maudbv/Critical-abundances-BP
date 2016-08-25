#### Geographical clustering/spatial autocorrelation of impacted plots

#############       Using GIS data
require(rgdal)

study_area <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_study_area")
geology <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_geology_map")
population <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_population_data")
plots <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_points")

coordinates(envplot) <- c("POINTX", "POINTY")
proj4string(envplot) <- proj4string(study_area)

envplot.grid <- points2grid(envplot, tolerance = 0.63, round=1)
grid <- SpatialGrid(envplot.grid, proj4string = proj4string(study_area))

# geographical distance matrix
library(raster)
spDst <-pointDistance(coordinates(envplot), lonlat = F, allpairs=FALSE)
rownames(spDst) <- colnames(spDst) <-rownames(coordinates(envplot))


# autocorrelation among impacted and non impacted sites

## home made null model of spatial mean nearest distance
spatial.mnnd <- function(sp.dist = spDst, X = var==2, nreps = 999, plot.histo = F) {

stopifnot(dim(sp.dist)[1] == length(X))
diag(sp.dist) <- NA
mnnd.obs <- mean(apply(sp.dist[X, X], 1, min, na.rm=T))

mnnd.rand = rep(NA, nreps)
for (i in 1:nreps){
  rv <- sample(X)
  mnnd.rand[i] <- mean(apply(sp.dist[rv, rv], 1, min, na.rm=T))
}

mnnd.rand <-  c(mnnd.obs,  mnnd.rand)
if (plot.histo) {
hist(mnnd.rand)
abline(v = mnnd.obs, col ="red")
}
p.mnnd <- (sum(mnnd.rand <mnnd.obs) + sum(mnnd.rand ==mnnd.obs)/2)/(nreps+1)


return(cbind(obs = mnnd.obs, mean.null = mean(mnnd.rand), sd.null = sd(mnnd.rand),P = p.mnnd))
}

### apply to target species threshold  ######

mnnd.th <- matrix(NA, nrow= 11, ncol = 8)
colnames(mnnd.th) <- paste(c("mnnd.obs","null.mean", "null.sd", "P<obs"),
               c(rep("below", 4), rep("above", 4)), sep = "_")
rownames(mnnd.th) <- impsp

for (i in 1:length(impsp)) {

sp=impsp[i]


community <- comm
community <- community[!is.na(community[,sp]) & community[,sp]>0,]  #only plots where the species is occurring

dst <- spDst[rownames(community),rownames(community)]   # select the corresponding distance submatrix

var <- community[,sp]  # extract the vector of abundances per plot
var[var < glmSRnat.overall$impact.spread [sp,"th.CI"] & var!=0] <- 1
var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 2


mnnd.th [i,1:4] <- spatial.mnnd(sp.dist = dst, X = var==1, nreps = 999)
mnnd.th [i,5:8] <- spatial.mnnd(sp.dist = dst, X = var==2, nreps = 999)
print(paste(i, sp))
}


################### ADDITIONAL spatial descriptions

#### apply to all abundance levels * trends of distinct groups of plots * ##########

tmp<- matrix(NA, nrow= 11, ncol = 6)
colnames(tmp) <- paste("c", 1:6, sep = "")
rownames(tmp) <- impsp
mnnd.trends <- list(obs = tmp, mean.null = tmp, sd.null = tmp, z = tmp, P = tmp)

for (i in 1:length(impsp)) {

  sp=impsp[i]

  community <- comm
  community <- community[!is.na(community[,sp]) & community[,sp]>0,]
  dst <- spDst[rownames(community),rownames(community)]

  var <- community[,sp]
  abun = as.numeric(names(table(var)[table(var)>5]))

  for (k in abun) {   ### loop on all the levels of abundance

  mnnd <- spatial.mnnd(sp.dist = dst, X = var==k, nreps = 999, plot.histo = F)

  mnnd.trends$obs[i,k] <- mnnd[,"obs"]
  mnnd.trends$mean.null[i,k] <- mnnd[,"mean.null"]
  mnnd.trends$sd.null[i,k] <- mnnd[,"sd.null"]
  mnnd.trends$z[i,k] <- (mnnd[,"obs"] - mnnd[,"mean.null"])/ mnnd[,"sd.null"]
  mnnd.trends$P[i,k] <- mnnd[,"P"]

  print(paste(i, sp, "class", k))

  }
}

# ### graph
#
# plot(mnnd.trends$z[1,], type= "b", ylim = c(-4,2), xlim=c(0,6), ann=F,pch = 21,
#      bg = c(NA, "black")[(mnnd.trends$P[1,]<=0.05) +1])
# text(1,mnnd.trends$z[1,1], pos = 2, label= impsp[1], cex = 0.5 )
# for( i in 2:11) {
#   par(new=T)
#   plot(mnnd.trends$z[i,], type= "b", ylim =c(-4,2), xlim=c(0,6), ann= F, axes=F,pch = 21,
#        bg = c(NA, "black")[(mnnd.trends$P[i,]<=0.05) +1])
#   text(1,mnnd.trends$z[i,1], pos = 2, label= impsp[i], cex = 0.5 )
# }
# abline(h=0, lty="dotted")
# mtext(2, text = "SES(MNND) of plots == abundance", line=2.5)
# mtext(1, text = "abundance class", line=2.5)

#### apply to all abundance levels *hierarchical groups ABOVE abundances* ##############

tmp<- matrix(NA, nrow= 11, ncol = 6)
colnames(tmp) <- paste("c", 1:6, sep = "")
rownames(tmp) <- impsp
mnnd.above <- list(obs = tmp, mean.null = tmp, sd.null = tmp, z = tmp,P = tmp)

for (i in 1:length(impsp)) {
  sp=impsp[i]

  community <- comm
  community <- community[!is.na(community[,sp]) & community[,sp]>0,]
  dst <- spDst[rownames(community),rownames(community)]

  var <- community[,sp]
  abun = as.numeric(names(table(var)[table(var)>5]))

  for (k in abun) {   ### loop on all the levels of abundance

    mnnd <- spatial.mnnd(sp.dist = dst, X = var>=k, nreps = 999, plot.histo = F)

    mnnd.above$obs[i,k] <- mnnd[,"obs"]
    mnnd.above$mean.null[i,k] <- mnnd[,"mean.null"]
    mnnd.above$sd.null[i,k] <- mnnd[,"sd.null"]
    mnnd.above$z[i,k] <- (mnnd[,"obs"] - mnnd[,"mean.null"])/ mnnd[,"sd.null"]
    mnnd.above$P[i,k] <- mnnd[,"P"]

    print(paste(i, sp, "class", k))

  }
}

#
# ### graph
# plot(mnnd.above$z[1,], type= "b", ylim = c(-5,2), xlim=c(1,6), ann=F, pch = 21,
#      bg = c(NA, "black")[(mnnd.above$P[1,]<=0.05) +1])
# text(2,mnnd.above$z[1,2], pos = 2, label= impsp[1], cex = 0.5 )
#
# for( i in 2:11) {
#   par(new=T)
#   plot(mnnd.above$z[i,], type= "b", ylim =c(-5,2), xlim=c(1,6), ann= F, axes=F, pch = 21,
#        bg = c(NA, "black")[(mnnd.above$P[i,]<=0.05) +1])
#   text(2,mnnd.above$z[i,2], pos = 2, label= impsp[i], cex = 0.5 )
# }
# abline(h=0, lty="dotted")
# mtext(2, text = "SES(MNND) of plots > abundance", line=2.5)
# mtext(1, text = "abundance class", line=2.5)



########## GRAPHICAL OUTPUT ##########

#  spatial clustering vs each abundance
#
# par(mfrow = c(3,4), mar=c(0,0,2,2), oma=c(6,6,2,3))
# ylim=c(-5,2)
# for (i in 1 : length(impsp)){
#
#   if (i ==4) plot.new()
#   sp <- impsp[i]
#   # identify significant classes in black
#   cols <- c(NA, "black") [ ( mnnd.above$P[i,] <=0.025) +1]
#   n <- as.numeric(M$n.obs[sp,])[2:6]
#   # create x axis for the 5 abundance classes
#   x = c(0:5)
#
#   # create y axis with the gamma richness Standardized effect size
#   y = mnnd.above$z[i,]
#   y [ (2:6)[ n < 5] ]=NA
#   y[1]=NA
#
#
#   # plot background
#   plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
#        pch = 21, col ="black",bg = cols)
#
#   axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
#   if (i %in% c(8,9,10,11)) text(y=-5.2, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)
#
#   if ( i %in% c(1,4,8)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
#   abline(h=0,lty="dotted")
#
#   #plot points and lines
#   par(new=T)
#   plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)
#
#   # Add species name
#   mtext(3, text=paste(species[sp, "Genus"]," ", species[sp, "Species"], sep="") ,
#         font = 3, outer= F,adj=0.9, cex=0.7, line=0.2, las = 1)
#   #   # Y axis label
#   #   if ( i %in% c(1,4,7)) {
#   #     mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
#   #   }
# }
#
# mtext(1, text=c("Abundance class"), adj=0.5, line=4, las = 1, outer=T)
# mtext(2, text=c("SES(MNND) of plots > abundance"), adj=0.5, line=3, las = 0, outer=T)
# #########  OTHER autocorrelation methods to be pursued #########
#
# ### mantel correlogram
# library(ecodist)
#
# # analyze the pattern of z across space
# space.d <- as.dist(spDst)
# SRnat.d <- dist(envplot$SRnat, "eucl")
# SRnat.mgram <- mgram(SRnat.d, space.d, nperm=100, nclass =10)
# plot(SRnat.mgram, main = "SRnat")
#
# x11()
# space.d <- as.dist(spDst)
# SRali.d <- dist(envplot$SRali, "eucl")
# SRali.mgram <- mgram(SRali.d, space.d, nperm=100, nclass =10)
# plot(SRali.mgram, main = "SRali")
#
# #### spatial autocor of abundances with distances for 11 impsp
# ab.mgram <- as.list(impsp)
# names(ab.mgram) <- impsp
# for (sp in impsp){
# community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
# ab <- community[,sp]
# plots<- rownames(community)
# space.d <- as.dist(spDst[plots,plots])
# ab.d <- dist(ab, "eucl")
# f <- mgram(ab.d, space.d, nperm=100, nclass =5)
# ab.mgram[[sp]] <- f
# }
#
# x11()
# par(mfrow= c(3,4))
# for (sp in impsp){
# f <- ab.mgram[[sp]]
# plot(f, main=sp)
# }
#
# ####### Beta dissimilarity in composition with distance overall   #########
#
#   community <- comm[which(rownames(comm) %in% realgrasslands),colnames(comm) %in% natives]
#   community <-community[rowSums(community) > 0,]
#   community[] <- as.numeric(community>0)
#
#   beta <- vegdist(community)
#   plots<- rownames(community)
#   space.d <- as.dist(spDst[plots,plots])
#   beta.d <- dist(beta, "eucl")
#   f <- mgram(beta.d, space.d, nperm=100, nclass =10)
#   plot(f)
#  abline(h=0)
#
#
#
#
# beta.mgram <-as.list(impsp)
# par(mfrow= c(3,4))
# for ( i in 1:length(impsp)) {
#
# sp <- impsp[i]
# community <- comm[which(rownames(comm) %in% realgrasslands & comm[,sp]>0 ),colnames(comm) %in% natives]
# var <- comm[rownames(community),sp]
# abun <-  unique(var)
#
# plot(range(spDst), c(-0.2, 0.2), type = "n", xlim =range(spDst), ylim= c(-0.2, 0.2))
# abline(h=0)
# beta.mgram[[sp]] <-  lapply(1:6, FUN = function(k) {
# community <- comm[which(rownames(comm) %in% realgrasslands & (comm[,sp]==k)),colnames(comm) %in% natives]
# # community <-community[rowSums(community) > 0,]
# f <- NA
# if (sum(rowSums(community) >0) >5) {
# beta <- vegdist(community, method = "euclid")
# plots<- rownames(community)
# space.d <- as.dist(spDst[plots,plots])
# beta.d <- dist(beta, "eucl")
# f <- mgram(beta.d, space.d, nperm=0, nclass =5)
# par(new = T)
# plot(f, ann = F, axes = F, xlim =range(spDst), ylim=c(-0.2, 0.2), main =sp,
#      col = paste("grey", seq(0,100,20),sep ="")[k] )
# }
# return(f)
# })
#
# }
