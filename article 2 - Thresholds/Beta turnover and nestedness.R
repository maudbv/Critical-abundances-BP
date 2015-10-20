
# Beta diversity partitioning
betasor<- function(spnames = impsp, com = comm, group = natives, nreps = 99) {

  betasor.list    <- list()

  for (i in 1:length(spnames)) {

    sp = spnames[i]
    if (spnames == "overall") {
      com$fake <- 1
      sp = "fake"}


    community <- com[which((rownames(com) %in% realgrasslands) & (com[,sp]>0) ),]
    community <- community[,colSums(community)>0]
    community<- community[,colnames(community) %in% group]
    community<- ceiling(community>0)

    var <- com[rownames(community),sp]
    var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
    var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1


    # beta partitioning
    nbeta.below <- (if(sum(var==0) == 0) {NULL}
                    else {
                      oecosimu(community[var==0,], nestedbetasor, method = "r1", nsimul = nreps, burnin = 0, thin = 1,
                            statistic = "statistic", alternative = c("less"))})

    nbeta.above <- (if(sum(var==1) == 0) {NULL}
                    else {
                      oecosimu(community[var==1,], nestedbetasor, method = "r1", nsimul = nreps, burnin = 0, thin = 1,
                            statistic = "statistic", alternative ="less")})


    betasor.list[[i]] <- cbind(below = ( if (!is.null(nbeta.below) ) { data.frame(obs = nbeta.below$oecosimu$statistic,
                                                  mean=nbeta.below$oecosimu$means,
                                                  z=nbeta.below$oecosimu$z,
                                                  pval = nbeta.below$oecosimu$pval)}
                               else {data.frame(obs = rep(NA,3),
                                                mean=rep(NA,3),
                                                z=rep(NA,3),
                                                pval = rep(NA,3)) }),
                               above = data.frame(obs = nbeta.above$oecosimu$statistic,
                                                  mean=nbeta.above$oecosimu$means,
                                                  z=nbeta.above$oecosimu$z,
                                                  pval = nbeta.above$oecosimu$pval))


    names(betasor.list)[i] <- sp

  }
  return(betasor.list)

}



# Beta diversity partitioning = my version with null model mixing above and below obs.
betasor<- function(spnames = impsp, com = comm, group = natives, nreps = 99) {

  betasor.list    <- list()

  for (i in 1:length(spnames)) {

    sp = spnames[i]

    # select plots of grassland where focal sp is present
    community <- com[which((rownames(com) %in% realgrasslands) & (com[,sp]>0) ),]
    community <- community[,colSums(community)>0] # remove species never present
    community<- community[,colnames(community) %in% group] # keep only native or alien species
    community<- ceiling(community>0) # make binary (presence / absence)
    community <- community[rowSums(community) >0,] # remove empty plots

    # group plots according to the level of abundance of the focal species
    var <- com[rownames(community),sp]
    var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0 # below critical level = 0
    var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1 # above or equal to critical level = 1

    # beta partitioning
    obs.below <- nestedbetasor(community[var==0,])
    obs.above <- nestedbetasor(community[var==1,])

    # permutation test
  n.below <-   n.above <- matrix(NA, ncol =3, nrow = nreps,
                                 dimnames = list(1:nreps, c("turnover", "nestedness","sorensen" )))
  for (j in 1:nreps) {
    rcommunity <- simulate(nullmodel(community, method ="r1"), nsim =1)
    random_plots <- sample(var, replace =FALSE)
    n.below[j,] <- nestedbetasor(community[ random_plots == 0 , ])
    n.above[j,] <- nestedbetasor(community[ random_plots == 1 , ])
  }

  n.below <- rbind(obs.below, n.below)
  n.above <- rbind(obs.above, n.above)

p.above <- (colSums(apply(n.above, 2, function(x) x<x[1]))
            + colSums(apply(n.above, 2, function(x) x==x[1]))/2) / nreps
p.below <- (colSums(apply(n.below, 2, function(x) x<x[1]))
            + colSums(apply(n.below, 2, function(x) x==x[1]))/2) / nreps

mean.above <- apply(n.above, 2, mean)
mean.below <- apply(n.below, 2, mean)

sd.above <- apply(n.above, 2, sd)
sd.below <- apply(n.below, 2, sd)

z.above <- apply(n.above, 2, function(x) (x[1] - mean(x))/sd(x))
z.below <- apply(n.below, 2, function(x) (x[1] - mean(x))/sd(x))

betasor.list[[i]] <-cbind( obs.below,mean.below,sd.below, z.below, p.below,
                           obs.above,mean.above,sd.above, z.above, p.above)
names(betasor.list)[i] <- sp
  }
  return(betasor.list)

}



# Function for Beta nestedness dissimilarity using BAselga et al 2010 :

betasor.dist<- function(com = community) {
  require(vegan)
  sor.dist <- betadiver(com, method = "w")   # (b+c) / (2a + b + c)
  sim.dist <- betadiver(com, method = "sim")   # pmin(b+c) / (a + pmin(b + c))
  nest.dist <- sor.dist - sim.dist
  betadist.list <- list(beta.sorensen =sor.dist, beta.simpson = sim.dist , beta.nested = nest.dist)
}


# Beta diversity partitioning = my version with null model mixing above and below obs.
beta.dissim.part<- function(sp="ACHMIL", com = comm, group = natives, nreps = 99) {

  betasor.list    <- list()

  for (i in 1:length(spnames)) {

    sp = spnames[i]

    # select plots of grassland where focal sp is present
    community <- com[which((rownames(com) %in% realgrasslands) & (com[,sp]>0) ),]
    community <- community[,colSums(community)>0] # remove species never present
    community<- community[,colnames(community) %in% group] # keep only native or alien species
    community<- ceiling(community>0) # make binary (presence / absence)
    community <- community[rowSums(community) >0,] # remove empty plots

    # group plots according to the level of abundance of the focal species
    var <- com[rownames(community),sp]
    var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0 # below critical level = 0
    var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1 # above or equal to critical level = 1

    # beta partitioning
    obs.below <- nestedbetasor(community[var==0,])
    obs.above <- nestedbetasor(community[var==1,])

    # permutation test
    n.below <-   n.above <- matrix(NA, ncol =3, nrow = nreps,
                                   dimnames = list(1:nreps, c("turnover", "nestedness","sorensen" )))
    for (j in 1:nreps) {
      rcommunity <- simulate(nullmodel(community, method ="r1"), nsim =1)
      random_plots <- sample(var, replace =FALSE)
      n.below[j,] <- nestedbetasor(community[ random_plots == 0 , ])
      n.above[j,] <- nestedbetasor(community[ random_plots == 1 , ])
    }

    n.below <- rbind(obs.below, n.below)
    n.above <- rbind(obs.above, n.above)

    p.above <- (colSums(apply(n.above, 2, function(x) x<x[1]))
                + colSums(apply(n.above, 2, function(x) x==x[1]))/2) / nreps
    p.below <- (colSums(apply(n.below, 2, function(x) x<x[1]))
                + colSums(apply(n.below, 2, function(x) x==x[1]))/2) / nreps

    mean.above <- apply(n.above, 2, mean)
    mean.below <- apply(n.below, 2, mean)

    sd.above <- apply(n.above, 2, sd)
    sd.below <- apply(n.below, 2, sd)

    z.above <- apply(n.above, 2, function(x) (x[1] - mean(x))/sd(x))
    z.below <- apply(n.below, 2, function(x) (x[1] - mean(x))/sd(x))

    betasor.list[[i]] <-cbind( obs.below,mean.below,sd.below, z.below, p.below,
                               obs.above,mean.above,sd.above, z.above, p.above)
    names(betasor.list)[i] <- sp
  }
  return(betasor.list)

}

