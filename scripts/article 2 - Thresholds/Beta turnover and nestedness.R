


# Beta diversity partitioning = my version with null model mixing above and below obs.
## Uses funcitons: nestedbetasor ( Baselga (2010)), from vegan, and simulate, from stats
betasor.multi<- function(spnames = impsp, com = comm, group = natives, null.comm = T, bootstrap = T, nreps = 99) {

  betasor.list    <- list()
require(vegan)
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
    obs.diff <- obs.above - obs.below
      
    # permutation test
    n.below <-   n.above <- matrix(NA, ncol =3, nrow = nreps,  dimnames = list(1:nreps, c("turnover", "nestedness","sorensen" )))
  for (j in 1:nreps) {
    rcommunity <- community
    rep.below <- which(var == 0)
    rep.above <- which(var == 1)
    
    if( null.comm == T) {
    # randomizes species associations within communities, and thus beta diversity. 
    # keeps SR and species freauency constant across landscape.
     rcommunity <- as.data.frame(simulate(nullmodel(community, method ="r1"), nsim =1)[])
    }
    
    if( bootstrap == T) {
     # + bootstrap plot above or below abundance (randomizes species richness above and below)
     rep.below <- sample(rep.below, replace =TRUE)
     rep.above <- sample(rep.above, replace =TRUE)
    }
    
    n.below[j,] <- nestedbetasor(rcommunity[ rep.below , ])
    n.above[j,] <- nestedbetasor(rcommunity[ rep.above , ])
  }
  n.below <- rbind(obs.below, n.below)
  n.above <- rbind(obs.above, n.above)
  n.diff <-   n.above  -   n.below 

## two sided p value:  
calc.pval <- function(x, n = nreps) (min(sum(x<x[1]), sum(x>x[1])) + sum(x==x[1])) /(n+1)

p.above <- apply(n.above, 2, calc.pval)
p.below <- apply(n.below, 2, calc.pval)
p.diff <- apply(n.diff, 2, calc.pval)

mean.above <- apply(n.above, 2, mean)
mean.below <- apply(n.below, 2, mean)
mean.diff <- apply(n.diff, 2, mean)

sd.above <- apply(n.above, 2, sd)
sd.below <- apply(n.below, 2, sd)
sd.diff <- apply(n.diff, 2, sd)

q2.5.above <- apply(n.above, 2, quantile, 0.025)
q2.5.below <- apply(n.below,  2, quantile, 0.025)
q2.5.diff <- apply(n.diff, 2, quantile, 0.025)

q97.5.above <- apply(n.above, 2, quantile, 0.975)
q97.5.below <- apply(n.below,  2, quantile, 0.975)
q97.5.diff <- apply(n.diff, 2, quantile, 0.975)

z.above <- apply(n.above, 2, function(x) (x[1] - mean(x))/sd(x))
z.below <- apply(n.below, 2, function(x) (x[1] - mean(x))/sd(x))
z.diff <- apply(n.diff, 2, function(x) (x[1] - mean(x))/sd(x))

betasor.list[[i]] <-cbind( obs.below, q2.5.below, q97.5.below, mean.below,sd.below, z.below, p.below,
                           obs.above, q2.5.above, q97.5.above, mean.above,sd.above, z.above, p.above,
                           obs.diff, q2.5.diff, q97.5.diff, mean.diff,sd.diff, z.diff, p.diff)
names(betasor.list)[i] <- sp
  }
  return(betasor.list)

}

# OTHER Function for Beta nestedness dissimilarity using BAselga et al 2010 :
betasor.dist<- function(com = community, group = natives, null.comm = T, bootstrap = T, nreps = 99) {
  
  require(vegan)
  
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
  
#calculate observed values
  sor.dist <- betadiver(community, method = "w")   # (b+c) / (2a + b + c)
  sim.dist <- betadiver(community, method = "sim")   # pmin(b+c) / (a + pmin(b + c))
  nest.dist <- sor.dist - sim.dist
  
  # Bootstrap resampling
  n.below <-   n.above <- matrix(NA, ncol =3, nrow = nreps,  dimnames = list(1:nreps, c("turnover", "nestedness","sorensen" )))
  for (j in 1:nreps) {
## TO DO ##
    }
    
    n.below[j,] <- nestedbetasor(rcommunity[ rep.below , ])
    n.above[j,] <- nestedbetasor(rcommunity[ rep.above , ])
  }
 
 n.below <- rbind(obs.below, n.below)
  n.above <- rbind(obs.above, n.above)
  n.diff <-   n.above  -   n.below 
  
  ## two sided p value:  
  calc.pval <- function(x, n = nreps) (min(sum(x<x[1]), sum(x>x[1])) + sum(x==x[1])) /(n+1)
  
  p.above <- apply(n.above, 2, calc.pval)
  p.below <- apply(n.below, 2, calc.pval)
  p.diff <- apply(n.diff, 2, calc.pval)
  
  mean.above <- apply(n.above, 2, mean)
  mean.below <- apply(n.below, 2, mean)
  mean.diff <- apply(n.diff, 2, mean)
  
  sd.above <- apply(n.above, 2, sd)
  sd.below <- apply(n.below, 2, sd)
  sd.diff <- apply(n.diff, 2, sd)
  
  q2.5.above <- apply(n.above, 2, quantile, 0.025)
  q2.5.below <- apply(n.below,  2, quantile, 0.025)
  q2.5.diff <- apply(n.diff, 2, quantile, 0.025)
  
  q97.5.above <- apply(n.above, 2, quantile, 0.975)
  q97.5.below <- apply(n.below,  2, quantile, 0.975)
  q97.5.diff <- apply(n.diff, 2, quantile, 0.975)
  
  z.above <- apply(n.above, 2, function(x) (x[1] - mean(x))/sd(x))
  z.below <- apply(n.below, 2, function(x) (x[1] - mean(x))/sd(x))
  z.diff <- apply(n.diff, 2, function(x) (x[1] - mean(x))/sd(x))
  
  betasor.list[[i]] <-cbind( obs.below, q2.5.below, q97.5.below, mean.below,sd.below, z.below, p.below,
                             obs.above, q2.5.above, q97.5.above, mean.above,sd.above, z.above, p.above,
                             obs.diff, q2.5.diff, q97.5.diff, mean.diff,sd.diff, z.diff, p.diff)
  names(betasor.list)[i] <- sp

  betadist.list <- list(beta.sorensen =sor.dist, beta.simpson = sim.dist , beta.nested = nest.dist)
}

# OTHER VERSION (not tested) Beta diversity partitioning using oecosimu
# betasor<- function(spnames = impsp, com = comm, group = natives, nreps = 99) {
#   
#   betasor.list    <- list()
#   
#   for (i in 1:length(spnames)) {
#     
#     sp = spnames[i]
#     if (spnames == "overall") {
#       com$fake <- 1
#       sp = "fake"}
#     
#     
#     community <- com[which((rownames(com) %in% realgrasslands) & (com[,sp]>0) ),]
#     community <- community[,colSums(community)>0]
#     community<- community[,colnames(community) %in% group]
#     community<- ceiling(community>0)
#     
#     var <- com[rownames(community),sp]
#     var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
#     var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
#     
#     
#     # beta partitioning
#     nbeta.below <- (if(sum(var==0) == 0) {NULL}
#                     else {
#                       oecosimu(community[var==0,], nestedbetasor, method = "r1", nsimul = nreps, burnin = 0, thin = 1,
#                                statistic = "statistic", alternative = c("less"))})
#     
#     nbeta.above <- (if(sum(var==1) == 0) {NULL}
#                     else {
#                       oecosimu(community[var==1,], nestedbetasor, method = "r1", nsimul = nreps, burnin = 0, thin = 1,
#                                statistic = "statistic", alternative ="less")})
#     
#     
#     betasor.list[[i]] <- cbind(below = ( if (!is.null(nbeta.below) ) { data.frame(obs = nbeta.below$oecosimu$statistic,
#                                                                                   mean=nbeta.below$oecosimu$means,
#                                                                                   z=nbeta.below$oecosimu$z,
#                                                                                   pval = nbeta.below$oecosimu$pval)}
#                                          else {data.frame(obs = rep(NA,3),
#                                                           mean=rep(NA,3),
#                                                           z=rep(NA,3),
#                                                           pval = rep(NA,3)) }),
#                                above = data.frame(obs = nbeta.above$oecosimu$statistic,
#                                                   mean=nbeta.above$oecosimu$means,
#                                                   z=nbeta.above$oecosimu$z,
#                                                   pval = nbeta.above$oecosimu$pval))
#     
#     
#     names(betasor.list)[i] <- sp
#     
#   }
#   return(betasor.list)
#   
# }
# 
