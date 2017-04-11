### Multivariate analysis of the invaded communities composition

library(vegan)
require(parallel)
require(doParallel)

# # NMDS all species in all grassland sites
community <- comm[which(rownames(comm) %in% realgrasslands),] 
community <- community[,colSums(community)>0]  ## 466 species

start <- cmdscale(vegdist(community, dist = "jaccard"),  k = 2)

clus <- makeCluster(6)
registerDoParallel(clus)
clusterEvalQ(clus, library(vegan))
nmds.complete <- metaMDS(community, distance = "jaccard", k = 2, trymax = 100, 
                         previous.best = start, parallel=6)
stopCluster(cl)

nmds.results <- nmds.complete
lim=c(min(nmds.results$points, na.rm=T), max(nmds.results$points, na.rm = T))
plot(nmds.results$points, col="darkgrey",ylim=lim, xlim=lim )
x <- nmds.results$species[impsp,1]
y <- nmds.results$species[impsp,2]
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7 )
mtext(3, text ="complete communities")

outliers <- identify(nmds.results$points)
outliers <- rownames(nmds.results$points)[outliers]

#remove outliers:
community <- community [!rownames(community)%in% outliers,]

start <- cmdscale(vegdist(community, dist = "jaccard"),  k = 2)

clus <- makeCluster(6)
registerDoParallel(clus)
clusterEvalQ(clus, library(vegan))
nmds.complete.noout <- metaMDS(community, distance = "jaccard", k = 2, trymax = 100, 
                         previous.best = start, parallel=6)
stopCluster(clus)


nmds.results <- nmds.complete.noout
lim=c(min(nmds.results$points, na.rm=T), max(nmds.results$points, na.rm = T))
plot(nmds.results$points, col="darkgrey",ylim=lim, xlim=lim )
x <- nmds.results$species[impsp,1]
y <- nmds.results$species[impsp,2]
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7 )
mtext(3, text ="complete communities")
# plot gourps below and above threholds

nmds.results <- 



sp="ACHMIL"
  
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,colSums(community)>0]
  community<- community[,colnames(community) %in% natives] 
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1



# nmds.results$species[impsp]




# NMDS of only target species co-occurence
# nmds.impsp <- list()
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
sp=impsp[i]
community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
community <- community[,names(community) %in% natives] 
community <- ceiling(community>0)


var <- comm[rownames(community),sp]
var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1

var2 <-as.character(var)

fact <- matrix(paste(matrix(var2, length(var), length(var)) ,
             matrix(var2, length(var), length(var), byrow=T), sep="-"),
       nrow = length(var), ncol= length(var))

D <- d.SM.SRnat[[i]]
tm <- lower.tri(fact)
plot(as.matrix(D)[tm] ~ as.factor(fact)[tm] )

# start <- cmdscale(D,  k = 2)
# nmds.impsp[[i]] <- metaMDS(D, k = 2, trymax = 10, previous.best = start)
#m <- nmds.impsp[[i]] 

# plot(m, display ="sites", type= "n")
# points( m$points[var==0,], col ="darkgrey")
# points( m$points[var==1,], pch =21, col ="red", bg ="red")
}

### Correlation matrix for the 11 species
community <- comm[which((rownames(comm) %in% realgrasslands) & (rowSums(comm[,impsp])>0) ),impsp] 

cor.impsp <- as.matrix(cor(community, method="spearman"))
diag(cor.impsp)<- NA
cor.impsp[lower.tri(cor.impsp)] <-NA
cor.impsp <-cor.impsp[,sort(impsp, decreasing = T)]

layout(matrix(c(rep(1,9), 2,2), nrow=11, ncol=1))
par(mar=c(3,5,2,2))
image(cor.impsp, col=colorRampPalette(c("blue","white","orangered"))(20),
      breaks=seq(-1, 1, 0.1), axes=F)
    axis(1, at =seq(0,1,0.1), label =rownames(cor.impsp), las=2,srt=30)
axis(2, at =seq(0,1,0.1), label =colnames(cor.impsp), las=1)

image.scale(cor.impsp, col=colorRampPalette(c("blue","white","orangered"))(20),
             breaks=seq(-1, 1, 0.1))     

cov.impsp <- cov(community)
diag(cov.impsp)<- NA
cov.impsp[lower.tri(cov.impsp)] <-NA
cov.impsp <-cov.impsp[,sort(impsp, decreasing = T)]

layout(matrix(c(rep(1,9), 2,2), nrow=11, ncol=1))
par(mar=c(3,5,2,2))
image(cov.impsp, col=colorRampPalette(c("blue","white","orangered"))(20),
      breaks=seq(-1, 1, 0.1), axes=F)
axis(1, at =seq(0,1,0.1), label =rownames(cov.impsp), las=2,srt=30)
axis(2, at =seq(0,1,0.1), label =colnames(cov.impsp), las=1)

image.scale(cov.impsp, col=colorRampPalette(c("blue","white","orangered"))(20),
            breaks=seq(-1, 1, 0.1))     


##PCA on 11 species ?


pca.impsp<- rda(community)


##### ORDINATION on presence - absence using RAUP CRICK

load("data for raupcrick.Rdata")

    memory.limit(4095)

    community <- comm[which((rownames(comm) %in% realgrasslands) ),]
    community <- community[,names(community) %in% natives] 
    
    dist.SR <- myraupcrick(community, nreps=999)
    dist <- raupcrick(community, nsimul=999) #it is a bit faster than mine but has pb with allocating memory
    gc()
    print(i)
    return(dist)
  })
)
stopCluster(cl)
save(dist.RC, file="saved Rdata/article 2 - threshold/dist.RC_Natives_99reps.Rdata")


## NMDS with Jaccard
community <- comm[which((rownames(comm) %in% realgrasslands) ),]
start <- cmdscale(vegdist(community, dist = "jaccard"),  k = 2)
nmds.SR <- metaMDS(community, distance = "jaccard", k = 3, trymax = 10, 
                               previous.best = start, parallel=6)

plot(nmds.SR)
outlier.sp <- identify(nmds.SR$species[,1:2], labels = row.names(nmds.SR$species))
# outlier.sp <- row.names(nmds.SR$species)[outlier.sp]

community <- comm[which((rownames(comm) %in% realgrasslands) ),- outlier.sp]
community <- community[rowSums(community)>0,]
community <- community[,colSums(community)>0]

start <- cmdscale(vegdist(community, dist = "jaccard"),  k = 2)
nmds.SR.clean<- metaMDS(community, distance = "jaccard", k = 3, trymax = 10, 
                               previous.best = start, parallel=6)
plot(nmds.SR.clean)
x <- nmds.SR.clean$species[impsp,1]*5
y <- nmds.SR.clean$species[impsp,2]*5
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.2,y +sign(y)*0.2, label = impsp, cex=0.8 )


community <- comm[which((rownames(comm) %in% realgrasslands) ),names(comm) %in% natives]
community <- community[,colSums(community)>0]
community <- community[rowSums(community)>0,]

start <- cmdscale(vegdist(community, dist = "raup"),  k = 2)
plot(start)
nmds.SRnat<- metaMDS(vegdist(community, dist = "raup"), k = 3, trymax = 10, 
                        previous.best = start, parallel=6)
plot(nmds.SRnat)


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


