### calculate community dissimilarities

save(impsp, comm, natives,aliens, realgrasslands, dissim.nm, myraupcrick, SMsim, file ="data for dissim nm.Rdata")

setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
library(parallel)
library(doParallel)

#########  SPECIES MATCHING raw values ##########
d.SM.SR<- list()

for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- ceiling(community>0)
  d.SM.SR[[i]] <- 1 - SMsim(community)
  }

d.SM.SRnat<- list()

for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  community <- ceiling(community>0)
  d.SM.SRnat[[i]] <- 1 - SMsim(community)
}

d.SM.SRali<- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% aliens] 
  community <- ceiling(community>0)
  d.SM.SR[[i]] <- 1 - SMsim(community)
}
save(d.SM.SR,d.SM.SRnat,d.SM.SRali, file = "saved Rdata/article 2 - threshold/d.SM.Rdata" )


#########  SPECIES MATCHING + NULL MODEL ###############

# Species matching + null model for ALL SPECIES
cl <- makeCluster(6)
system.time(
  dist.SM.SR <- parLapplyLB( cl = cl, X = 1:11, fun = function(i) {
#     setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R")
    load(file ="data for dissim nm.Rdata")
    library(vegan)
    memory.limit(4095)
    sp=impsp[i]
    # select plots of grassland where the target is present
    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
    community <- community[,which(colSums(community) >0) ]  # remove sp which never cooccur with sp
    
    dist <- dissim.nm(community, nullmodel = T, nreps =999, index = "sm")
    gc()
    print(i)
    return(dist)
  })
)
stopCluster(cl)
save(dist.SM.SR, file="saved Rdata/article 2 - threshold/dist.SM_allspecies.999reps.Rdata")


# Species matching + null model for NATIVES
cl <- makeCluster(6)
system.time(
  dist.SM.SRnat <- parLapplyLB( cl = cl, X = 1:11, fun = function(i) {
    #     setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R")
    load(file ="data for dissim nm.Rdata")
    library(vegan)
    memory.limit(4095)
    sp=impsp[i]
    # select plots of grassland where the target is present
    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
    community <- community[,names(community) %in% natives]  # keep only natives
    community <- community[,which(colSums(community) >0) ]  # remove natives which never cooccur with sp
    
    dist <- dissim.nm(community, nullmodel = T, nreps =499, index = "sm")
    gc()
    print(i)
    return(dist)
  })
)
stopCluster(cl)
save(dist.SM.SRnat, file="saved Rdata/article 2 - threshold/dist.SM_Natives_499reps.Rdata")



# Species matching + null model for ALIENS
cl <- makeCluster(6)
system.time(
  dist.SM.SRali <- parLapplyLB( cl = cl, X = 1:11, fun = function(i) {
    setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R")
    load(file ="data for dissim nm.Rdata")
    library(vegan)
    memory.limit(4095)
    sp=impsp[i]
    # select plots of grassland where the target is present
    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
    community <- community[,names(community) %in% aliens]  # keep only aliens
    community <- community[,which(colSums(community) >0) ]  # remove aliens which never cooccur with sp
    
    dist <- dissim.nm(community, nullmodel = T, nreps =499, index = "sm")
    gc()
    print(i)
    return(dist)
  })
)
stopCluster(cl)
save(dist.SM.SRali, file="saved Rdata/article 2 - threshold/dist.SM_aliens_499reps.Rdata")

## save all results
save(dist.SM.SR, dist.SM.SRnat,dist.SM.SRali, file="saved Rdata/article 2 - threshold/dist.SM_999reps.Rdata")
