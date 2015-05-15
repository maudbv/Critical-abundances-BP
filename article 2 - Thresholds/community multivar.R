### Multivariate analysis of the invaded communities composition

library(vegan)
require(parallel)
require(doParallel)

D <- glmSRnat.overall$impact.spread
comm[is.na(comm)] <- 0

# Species which show a threshold
impsp <- rownames(D[which(!is.na(D$th.CI) & (rownames(D) %in% aliens)),])



# # NMDS all species in all grassland sites
# community <- comm[which(rownames(comm) %in% realgrasslands),] 
# community <- community[,-which(colSums(community) == 0)]  ## 466 species
# 
# start <- cmdscale(vegdist(community, dist = "bray"),  k = 3)
# 
# clus <- makeCluster(6)
# registerDoParallel(clus)
# clusterEvalQ(clus, library(vegan))
# nmds.complete <- metaMDS(community, distance = "bray", k = 3, trymax = 100, 
#                          previous.best = start, parallel=6)
# stopCluster(cl)
# 
# nmds.results <- nmds.complete
# lim=c(min(nmds.results$points, na.rm=T), max(nmds.results$points, na.rm = T))
# plot(nmds.results$points, col="darkgrey",ylim=lim, xlim=lim )
# x <- nmds.results$species[impsp,1]
# y <- nmds.results$species[impsp,2]
# arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
# text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7 )
# mtext(3, text ="complete communities")
# 
# nmds.results$species[impsp]

# NMDS all species in sites with at least one target species
community <- comm[which((rownames(comm) %in% realgrasslands) & (rowSums(comm[,impsp])>0) ),] 
community <- community[,-which(colSums(community) == 0)] 

start <- cmdscale(vegdist(community, dist = "bray"),  k = 2)

clus <- makeCluster(6)
registerDoParallel(clus)
clusterEvalQ(clus, library(vegan))
nmds.complete.impsp <- metaMDS(community, distance = "bray", k = 3, trymax = 10, 
                               previous.best = start, parallel=6)
stopCluster(cl)

nmds.results <- nmds.complete.impsp
lim=c(min(nmds.results$points, na.rm=T), max(nmds.results$points, na.rm = T))
plot(nmds.results$points, col="grey",ylim=lim, xlim=lim , pch=20)
x <- nmds.results$species[impsp,1]*2
y <- nmds.results$species[impsp,2]*2
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.2,y +sign(y)*0.2, label = impsp, cex=0.8 )



# DECORANA all species in sites with at least one target species
community <- comm[which((rownames(comm) %in% realgrasslands) & (rowSums(comm[,impsp])>0) ),] 

decorana.complete.impsp <- decorana(community, iweigh=1, iresc=4, ira=0,
                                    mk=26, short=0)


quartz()
opar <- par(mfrow=c(1,2)) 
plot(decorana.complete.impsp , display="sites", col="grey", xlim=c(-3,3), ylim=c(-3,3))
x <- scores(decorana.complete.impsp , dis = "species")[impsp,1]
y <- scores(decorana.complete.impsp , dis = "species")[impsp,2]
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7 )
mtext(3, text ="DCA")

nmds.results <- nmds.complete.impsp
nmds.results$points[,1] <- - nmds.results$points[,1]
nmds.results$species[,1] <- - nmds.results$species[,1]


lim=c(-2, 2)
plot(nmds.results$points, col="darkgrey", xlim=lim, ylim=lim)
x <-  nmds.results$species[impsp,1]*2
y <- nmds.results$species[impsp,2]*2
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7 )
mtext(3, text ="NMDS")


##  EXAMPLE
### the detrending rationale:
gaussresp <- function(x,u) exp(-(x-u)^2/2)
x <- seq(0,6,length=15) ## The gradient
u <- seq(-2,8,len=23)   ## The optima
pack <- outer(x,u,gaussresp)
noisy <- (0.5 + runif(length(pack)))*pack
matplot(x, noisy, type="l", main="Species noisying")
opar <- par(mfrow=c(2,2))
plot(scores(cca(noisy ~ x), dis="sites"), asp=1, type="b", main="CCA")
plot(scores(capscale(noisy ~ x), dis="sites"), asp=1, type="b", main="CAPSCALE")
plot(scores(metaMDS(noisy, trymax=10)), asp=1, type="b", main="NMDS")
plot(scores(decorana(noisy)), asp=1, type="b", main="DCA")

# NMDS of only target species co-occurence
community <- comm[which((rownames(comm) %in% realgrasslands) & (rowSums(comm[,impsp])>0) ),impsp] 

start <- cmdscale(vegdist(community, dist = "bray"),  k = 2)
nmds.impsp <- metaMDS(community, distance = "bray", k = 2, trymax = 10, previous.best = start)

nmds.results <- nmds.impsp
lim=c(min(nmds.results$points, na.rm=T)+0.2, max(nmds.results$points, na.rm = T)+0.2)
plot(nmds.results$points, col="darkgrey",ylim=lim, xlim=lim )
x <- nmds.results$species[impsp,1]*2
y <- nmds.results$species[impsp,2]*2
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7 )
mtext(3, text ="only 11 aliens")


### heatmap of species distribution

hmap.impsp <- heatmap(apply(community, 2, as.numeric), distfun = vegdist,  method = "gower",
                      col = coltreat, labRow="")

cluster.impsp <-hclust(vegdist(t(apply(community, 2, as.numeric)), method = "gower"))




#### dbRDA per species abundance + ADONIS + BETADISPER
coltreat <-  colorRampPalette(c("palegoldenrod", "firebrick"))(7) 
coltreat.tr <- paste(coltreat, "50", sep="")


# dbrda.species<- list()
# adonis.species <- list()
# betadisper.species <- list()
# 
# for (i in 1:length(impsp)) {
#   sp=impsp[i]
#   community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
#   community <- community[,-which(colSums(community) == 0)] 
#   var <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),sp]
#   
#   dbrda.species[[i]]  <- capscale(community ~ var,  distance = "bray")  
#   adonis.species[[i]] <- adonis(community ~ var, method = "bray")
#   betadisper.species [[i]]<- betadisper(vegdist(community, method = "bray"), group = var,type ="median")
#   print(paste(i, ":", sp))
# }
# save(dbrda.species,adonis.species,betadisper.species, file = "saved Rdata/article 2 - threshold/dbRDA_bray.Rdata" )

#ordination plot
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.species[[i]]
  sp = impsp[i]
  
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  treat =  comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),sp]
  plot(dbrda.results, display = "sites", type="points", axes=F, ann=F) 
  
  for (j in sort(unique(treat))) {
    if (length(grep(j,treat)) >3) {
      #       huls <- ordihull(scores(dbrda.results, display = "sites")[grep(j,treat),],
      #        groups=treat[ treat == j],draw="polygon",
      #                        col=coltreat.tr[j+1],label=F)
      #    
    ordiellipse(scores(dbrda.results, display = "sites")[grep(j,treat),],
                  groups=treat[ treat == j],draw="polygon",display="sites",
                  kind = c("sd"), conf=0.95,
                  col=coltreat.tr[j+1], label=F, border="grey30")
      print(j)
      
    }
  }
  
  #   ordihull(dbrda.results,groups=treat,draw="polygon", label =T)
  f<- adonis.species [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.species [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}


plot(-1:8, -1:8, type ="n", ann=F, bty="n", axes=F)
for (i in 1:6) {
polygon(x = c(i,i,i+1,i+1),   y =c(4,5,5, 4), col =coltreat[i], border =F)
}
text(1.5:6.6,rep(3.5,6),
     label = 1:6, cex = 0.7)

text(1.5:6.6,rep(3,6), srt=45, adj=1,
    label = c("Occasional" ,"Common-Occasional",  "Common",
              "Abundant-Common", "Abundant","Dominant"), cex = 0.7)





#### dbRDA per ABOVE or below threshold
coltreat <-  colorRampPalette(c("palegoldenrod", "firebrick"))(7) 
coltreat.tr <- paste(coltreat, "50", sep="")


adonis.species.th <- list()
betadisper.species.th <- list()

for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,-which(colSums(community) == 0)] 
  var <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
   adonis.species.th[[i]] <- adonis(community ~ var, method = "bray")
  betadisper.species.th[[i]]<- betadisper(vegdist(community, method = "bray"), group = var,type ="median")
  print(paste(i, ":", sp))
}
save(dbrda.species,adonis.species,betadisper.species, file = "saved Rdata/article 2 - threshold/dbRDA_bray.Rdata" )

#ordination plot
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.species[[i]]
  sp = impsp[i]
  
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  var <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  treat <- c("below", "above") [var+1]
  plot(dbrda.results, display = "sites", type="none", pch ="", axes=F, ann=F) 
  
ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
             groups=treat[ treat == "below"],draw="polygon",
                             col="white",label=F)

ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
         groups=treat[ treat == "above"],draw="polygon",
         col="grey30",label=F)

#       ordiellipse(scores(dbrda.results, display = "sites")[grep(j,treat),],
#                   groups=treat[ treat == j],draw="polygon",display="sites",
#                   kind = c("sd"), conf=0.95,
#                   col=coltreat.tr[j+1], label=F, border="grey30")
       
  ordihull(dbrda.results,groups=treat,draw="polygon", label =T)
  f<- adonis.species.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.species.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}



#### dbRDA per ABOVE or below threshold for NATIVE RICHNESS
coltreat <-  colorRampPalette(c("palegoldenrod", "firebrick"))(7) 
coltreat.tr <- paste(coltreat, "50", sep="")

#
dbrda.th<- list()
adonis.th <- list()
betadisper.th <- list()

for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  if (sum(colSums(community) == 0)>0)   community <- community[,-which(colSums(community) == 0)] 
 if (sum(rowSums(community) == 0)>0) community <- community[-which(rowSums(community) == 0),] 
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1

  
  dbrda.th[[i]]  <- capscale(community ~ var,  distance = "bray")  
  adonis.th[[i]] <- adonis(community ~ var, method = "bray")
  betadisper.th [[i]]<- betadisper(vegdist(community, method = "bray"), group = var,type ="median")
  print(paste(i, ":", sp))
}
# save(dbrda.species,adonis.species,betadisper.species, file = "saved Rdata/article 2 - threshold/dbRDA_bray.Rdata" )

#ordination plot
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  if (sum(colSums(community) == 0)>0)   community <- community[,-which(colSums(community) == 0)] 
  if (sum(rowSums(community) == 0)>0) community <- community[-which(rowSums(community) == 0),] 
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  plot(dbrda.results, display = "sites", type="points", pch ="", axes=F, ann=F) 
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #       ordiellipse(scores(dbrda.results, display = "sites")[grep(j,treat),],
  #                   groups=treat[ treat == j],draw="polygon",display="sites",
  #                   kind = c("sd"), conf=0.95,
  #                   col=coltreat.tr[j+1], label=F, border="grey30")
  
  ordihull(dbrda.results,groups=treat,draw="polygon", label =T)
  f<- adonis.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}


#### raup crick distance  ?? => not working for now

distances <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,-which(colSums(community) == 0)] 
  
  var <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),sp]
  D<- raupcrick(community )
  distances[[i]] <- D
  print(paste(i, ":", sp))
}



