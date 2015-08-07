### Variations in Gamma diversity at each abundance for all species (or at least impsp to start with)
#   Input:
#         com: community matrix with rows = sites and columns = species codes
#         sp.names: vector of species codes for target species
#         group: vector of species codes (character strings) encompassing the component of diversity under study
#         nreps : integer. number of random sampling
#
# null.model: one of these :
#       permute.all : permute all sites randomly
#       permute.rare: permutes only the sites of a given abundance with sites where the species is rare
#       permute.below: permute only sites of a given abundance with all sites where abudnance is inferior
#
# output : a list
#         $obs : gamma trends for each species
#         $P: P value compared to null model for each species (quantile of null distribution)
#         $mean.null: mean null gamma value for each abundance class
#         $sd.null : standard deviation of null distribution
#         $q975: 97.5th quantile of null distribution
#         $q025: 2.5th quantile of null distribution

## REwrite of div part funciton with all abundance threshold explorations
alpha.trend <- function(spnames = impsp, group = natives, null.model = "permute.rare", nreps = 99) {

  # initiating output
  tmp<- as.data.frame(matrix(NA, nrow = length(spnames), ncol = 6,
                             dimnames= list(spnames,paste("c", 1: 6 , sep=""))))
  alpha.list     <- list( obs= tmp , P.alpha = tmp,mean.null= tmp, sd.null = tmp, q975= tmp , q025 = tmp )
  rm(tmp)

  # looping on each of the target species
  for (sp in spnames) {

    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
    community <- community[,colSums(community)>0]
    community<- community[,colnames(community) %in% group]
    community<- ceiling(community>0)

    var <- comm[rownames(community),sp]
    abun <- table(var)
    abun <- as.numeric(names(abun[abun>=5]))

    freq <- colSums(community)/dim(community)[1]
    rare.freq <- colSums(community[var==0,])/sum(var==0)

    alpha.obs <- as.data.frame(matrix(NA, nrow = 6, ncol = 6,
                                      dimnames= list(
                                        paste("c", 1: 6 ),
                                        c("nb.plot", "gamma", "alpha", "beta.add", "beta.mult","beta.prop")
                                      )))

    for (k in abun) {
      alpha.obs[k,] <-  c(n <- sum(rowSums(community[var==k,])>0),
                          gam <- sum(colSums(community[var==k,])>0),
                          alpha <- mean(rowSums(community[var==k,]), na.rm=T),
                          beta.add <- gam - alpha,
                          beta.mult <- gam/alpha,
                          beta.prop <- (gam - alpha)/gam
      )
    }

    # build random expectatinos
    r.alphas<- matrix(NA, nrow= nreps, ncol = 6)

    for (j in 1:nreps) {

      # null model A : simple permutation of all plots
      if( null.model == "permute.all") {
        rvar <- sample(var)
      }

      ### calculate alpha diversity at each abundance
      for (k in abun) {
        if( null.model == "permute.rare") {
          rvar <- var
          rvar[rvar%in% c(1, k)] <- sample(rvar[rvar%in% c(1, k)])
        }

        if( null.model == "permute.below") {
          rvar <- var
          rvar[rvar<= k] <- sample(rvar[rvar<= k])
        }

        r.alphas[j,k] <-  mean(rowSums(community[rvar==k,]), na.rm=T)

      }

      print(paste(sp, ":" ,j))

    }
    r.alphas <-  rbind( alpha.obs$alpha,  r.alphas)

    alpha.list$obs[sp,] <- alpha.obs$alpha
    alpha.list$P.alpha[sp,] <- apply(r.alphas, 2, function(x) (sum(x<x[1]) + sum(x == x[1])/2)/(nreps+1))
    alpha.list$mean.null[sp,] <-  apply(r.alphas, 2, function(x) mean(x, na.rm=T))
    alpha.list$sd.null[sp,] <-  apply(r.alphas, 2, function(x) sd(x, na.rm=T))
    alpha.list$q975[sp,] <-  apply(r.alphas, 2, function(x) quantile(x, probs = 0.975, na.rm=T))
    alpha.list$q025[sp,] <-  apply(r.alphas, 2, function(x) quantile(x, probs = 0.025, na.rm=T))


  }
  return(alpha.list)
}


## REwrite of div part funciton with all abundance threshold explorations
gamma.trend <- function(spnames = impsp, group = natives, null.model = "permute.rare", nreps = 99) {

# initiating output
  tmp<- as.data.frame(matrix(NA, nrow = length(spnames), ncol = 6,
                             dimnames= list(spnames,paste("c", 1: 6 , sep=""))))
  gamma.list     <- list( obs= tmp , alpha = tmp,
                          P.gamma = tmp,mean.null= tmp, sd.null = tmp, q975= tmp , q025 = tmp )
  rm(tmp)

# looping on each of the target species
  for (sp in spnames) {

    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
    community <- community[,colSums(community)>0]
    community<- community[,colnames(community) %in% group]
    community<- ceiling(community>0)

    var <- comm[rownames(community),sp]
    abun <- table(var)
    abun <- as.numeric(names(abun[abun>=5]))

    freq <- colSums(community)/dim(community)[1]
    rare.freq <- colSums(community[var==0,])/sum(var==0)

gamma.obs <- as.data.frame(matrix(NA, nrow = 6, ncol = 6,
                           dimnames= list(
                             paste("c", 1: 6 ),
                             c("nb.plot", "gamma", "alpha", "beta.add", "beta.mult","beta.prop")
                           )))

    for (k in abun) {
      gamma.obs[k,] <-  c(n <- sum(rowSums(community[var==k,])>0),
             gam <- sum(colSums(community[var==k,])>0),
             alpha <- mean(rowSums(community[var==k,]), na.rm=T),
             beta.add <- gam - alpha,
             beta.mult <- gam/alpha,
             beta.prop <- (gam - alpha)/gam
          )
    }

    # build random expectatinos
    r.gammas<- matrix(NA, nrow= nreps, ncol = 6)

    for (j in 1:nreps) {

      # null model A : simple permutation of all plots
      if( null.model == "permute.all") {
        rvar <- sample(var)
      }

      ### calculate gamma diversity at each abundance
     for (k in abun) {
              if( null.model == "permute.rare") {
                rvar <- var
                rvar[rvar%in% c(1, k)] <- sample(rvar[rvar%in% c(1, k)])
              }

              if( null.model == "permute.below") {
                rvar <- var
                rvar[rvar<= k] <- sample(rvar[rvar<= k])
              }

          r.gammas[j,k] <-  sum(colSums(community[rvar==k,])>0)
        }

     print(paste(sp, ":" ,j))

    }
   r.gammas <-  rbind( gamma.obs$gamma,  r.gammas)

   gamma.list$obs[sp,] <- gamma.obs$gamma
   gamma.list$alpha[sp,] <- gamma.obs$alpha
   gamma.list$P.gamma[sp,] <- apply(r.gammas, 2, function(x) (sum(x<x[1]) + sum(x == x[1])/2)/(nreps+1))
   gamma.list$mean.null[sp,] <-  apply(r.gammas, 2, function(x) mean(x, na.rm=T))
   gamma.list$sd.null[sp,] <-  apply(r.gammas, 2, function(x) sd(x, na.rm=T))
   gamma.list$q975[sp,] <-  apply(r.gammas, 2, function(x) quantile(x, probs = 0.975, na.rm=T))
   gamma.list$q025[sp,] <-  apply(r.gammas, 2, function(x) quantile(x, probs = 0.025, na.rm=T))


}
  return(gamma.list)
}


## REwrite of div part funciton with all abundance threshold explorations
beta.trend <- function(spnames = impsp, group = natives, null.model = "permute.rare", nreps = 99) {

  # initiating output
  tmp<- as.data.frame(matrix(NA, nrow = length(spnames), ncol = 6,
                             dimnames= list(spnames,paste("c", 1: 6 , sep=""))))
  beta.list     <- list( obs= tmp , P.beta = tmp, mean.null= tmp, sd.null = tmp, q975= tmp , q025 = tmp )
  rm(tmp)

  # looping on each of the target species
  for (sp in spnames) {

    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
    community <- community[,colSums(community)>0]
    community<- community[,colnames(community) %in% group]
    community<- ceiling(community>0)

    var <- comm[rownames(community),sp]

    freq <- colSums(community)/dim(community)[1]
    rare.freq <- colSums(community[var==0,])/sum(var==0)

    divpart.obs <- c(n <- sum(rowSums(community[var==1,])>0),
                   gam <- sum(colSums(community[var==1,])>0),
                   alpha <- mean(rowSums(community[var==1,]), na.rm=T),
                   beta.add <- gam - alpha,
                   beta.mult <- gam/alpha,
                   beta.prop <- (gam - alpha)/gam
    )

    for (k in 2: 6) {

      g <- ( if(sum(var==k)>1) {
        c(n <- sum(rowSums(community[var==k,])>0),
          gam <- sum(colSums(community[var==k,])>0),
          alpha <- mean(rowSums(community[var==k,]), na.rm=T),
          beta.add <- gam - alpha,
          beta.mult <- gam/alpha,
          beta.prop <- (gam - alpha)/gam
                )
      }
      else {
        c(n <- sum(sum(community[var==k,])>0),
          gam <- sum(community[var==k,]>0),
          alpha <- sum(community[var==k,], na.rm=T),
          beta.add <- 0,
          beta.mult <- 0,
          beta.prop <-0
        )
      })
      divpart.obs<- rbind(divpart.obs, g)
    }
    divpart.obs<- as.data.frame(divpart.obs)
    rownames(divpart.obs) <- paste("c", 1: 6 )
    colnames(divpart.obs) <- c("nb.plot", "gamma", "alpha", "beta.add", "beta.mult","beta.prop")


    # build random expectatinos
    r.betas<- matrix(NA, nrow= nreps, ncol = 6)

    for (j in 1:nreps) {

      # null model A : simple permutation of all plots
      if( null.model == "permute.all") {
        rvar <- sample(var)
      }

      ### calculate gamma diversity at each abundance
      r.betas[j,] <- sapply(1:6, FUN = function(k){

        if( null.model == "permute.rare") {
          rvar <- var
          rvar[rvar%in% c(1, k)] <- sample(rvar[rvar%in% c(1, k)])
        }

        if( null.model == "permute.below") {
          rvar <- var
          rvar[rvar<= k] <- sample(rvar[rvar<= k])
        }

        if (sum(rvar==k)>1) {
          # beta.prop
           b <- (sum(colSums(community[rvar==k,])>0) - mean(rowSums(community[var==k,]), na.rm=T))/sum(colSums(community[rvar==k,])>0)

#           #beta.miult
#           b <- sum(colSums(community[rvar==k,])>0)/ mean(rowSums(community[var==k,]), na.rm=T)

        }
        if (sum(rvar==k)<=1) b <-0
        return(b)
      })

      print(paste(sp, ":" ,j))
    }

    r.betas <-  rbind( divpart.obs$beta.prop,  r.betas)

    beta.list$obs[sp,] <- divpart.obs$beta.prop
    beta.list$P.beta[sp,] <- apply(r.betas, 2, function(x) (sum(x<x[1]) + sum(x == x[1])/2)/(nreps+1))
    beta.list$mean.null[sp,] <-  apply(r.betas, 2, function(x) mean(x, na.rm=T))
    beta.list$sd.null[sp,] <-  apply(r.betas, 2, function(x) sd(x, na.rm=T))
    beta.list$q975[sp,] <-  apply(r.betas, 2, function(x) quantile(x, probs = 0.975, na.rm=T))
    beta.list$q025[sp,] <-  apply(r.betas, 2, function(x) quantile(x, probs = 0.025, na.rm=T))
  }

  return(beta.list)
}

