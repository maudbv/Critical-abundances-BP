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
alpha.trend <-
  function(spnames = impsp, group = natives, null.model = "permute.rare", nreps = 99) {
    # initiating output
    tmp <- as.data.frame(matrix(
      NA, nrow = length(spnames), ncol = 6,
      dimnames = list(spnames,paste("c", 1:6 , sep =
                                      ""))
    ))
    alpha.list     <-
      list(
        obs = tmp , P.alpha = tmp,mean.null = tmp, sd.null = tmp, q975 = tmp , q025 = tmp
      )
    rm(tmp)

    # looping on each of the target species
    for (sp in spnames) {
      community <-
        comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp] > 0)),]
      community <- community[,colSums(community) > 0]
      community <- community[,colnames(community) %in% group]
      community <- ceiling(community > 0)

      var <- comm[rownames(community),sp]
      abun <- table(var)
      abun <- as.numeric(names(abun[abun >= 5]))

      freq <- colSums(community) / dim(community)[1]
      rare.freq <- colSums(community[var == 0,]) / sum(var == 0)

      alpha.obs <- as.data.frame(matrix(
        NA, nrow = 6, ncol = 6,
        dimnames = list(
          paste("c", 1:6),
          c(
            "nb.plot", "gamma", "alpha", "beta.add", "beta.mult","beta.prop"
          )
        )
      ))

      for (k in abun) {
        alpha.obs[k,] <-  c(
          n <- sum(rowSums(community[var == k,]) > 0),
          gam <- sum(colSums(community[var == k,]) > 0),
          alpha <-
            mean(rowSums(community[var == k,]), na.rm = T),
          beta.add <- gam - alpha,
          beta.mult <- gam / alpha,
          beta.prop <- (gam - alpha) / gam
        )
      }

      # build random expectatinos
      r.alphas <- matrix(NA, nrow = nreps, ncol = 6)

      for (j in 1:nreps) {
        # null model A : simple permutation of all plots
        if (null.model == "permute.all") {
          rvar <- sample(var)
        }


        ### calculate alpha diversity at each abundance
        for (k in abun) {
          if (null.model == "permute.rare") {
            rvar <- var
            rvar[rvar %in% c(1, k)] <- sample(rvar[rvar %in% c(1, k)])
          }

          if (null.model == "permute.below") {
            rvar <- var
            rvar[rvar <= k] <- sample(rvar[rvar <= k])
          }

          r.alphas[j,k] <-
            mean(rowSums(community[rvar == k,]), na.rm = T)

        }

        print(paste(sp, ":" ,j))

      }
      r.alphas <-  rbind(alpha.obs$alpha,  r.alphas)

      alpha.list$obs[sp,] <- alpha.obs$alpha
      alpha.list$P.alpha[sp,] <-
        apply(r.alphas, 2, function(x)
          (sum(x < x[1]) + sum(x == x[1]) / 2) / (nreps + 1))
      alpha.list$mean.null[sp,] <-
        apply(r.alphas, 2, function(x)
          mean(x, na.rm = T))
      alpha.list$sd.null[sp,] <-
        apply(r.alphas, 2, function(x)
          sd(x, na.rm = T))
      alpha.list$q975[sp,] <-
        apply(r.alphas, 2, function(x)
          quantile(x, probs = 0.975, na.rm = T))
      alpha.list$q025[sp,] <-
        apply(r.alphas, 2, function(x)
          quantile(x, probs = 0.025, na.rm = T))


    }
    return(alpha.list)
  }


## REwrite of div part funciton with all abundance threshold explorations
gamma.trend <-
  function(spnames = impsp, group = natives, bootstrapped = TRUE, null.model = "permute.rare",breps= 99, nreps = 99) {

     # initiating output
    tmp <- as.data.frame(matrix(
      NA, nrow = length(spnames), ncol = 6,
      dimnames = list(spnames,paste("c", 1:6 , sep =
                                      ""))
    ))

    gamma.list     <- list(
      obs = tmp , obs.deltagamma = tmp ,alpha = tmp,
      P.gamma = tmp,mean.null = tmp, sd.null = tmp, q975.null = tmp , q025.null = tmp,
      P.deltagamma = tmp,mean.null.deltagamma = tmp, sd.null.deltagamma = tmp, q975.null.deltagamma = tmp , q025.null.deltagamma = tmp,
      mean.boot = tmp, sd.boot = tmp, q975.boot = tmp , q025.boot = tmp
    )
    rm(tmp)

    # looping on each of the target species
    for (sp in spnames) {
      community <-comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp] > 0)),]
      community <- community[,colSums(community) > 0]
      community <- community[,colnames(community) %in% group]
      community <- ceiling(community > 0)

      var <- comm[rownames(community),sp]
      abun <- table(var)
      abun <- as.numeric(names(abun[abun >= 5]))
      freq <- colSums(community) / dim(community)[1]
      freq.rare <- colSums(community[var==1,]) / dim(community[var==1,])[1]

      gamma.obs <- as.data.frame(matrix(
        NA, nrow = 6, ncol = 6,
        dimnames = list(
          paste("c", 1:6),
          c(
            "nb.plot", "gamma", "alpha", "beta.add", "beta.mult","beta.prop"
          )
        )
      ))

      for (k in abun) {
        gamma.obs[k,] <-  c(
          n <- sum(rowSums(community[var == k,]) > 0),
          gam <- sum(colSums(community[var == k,]) > 0),
          alpha <- mean(rowSums(community[var == k,]), na.rm = T),
          beta.add <- gam - alpha,
          beta.mult <- gam / alpha,
          beta.prop <- (gam - alpha) / gam
        )
      }
      gamma.obs $deltagamma <-gamma.obs$gamma - gamma.obs$gamma[1]

      # BOOTSTRAPPING delta gamma
      if (bootstrapped) {
        # build random expectatinos
        b.deltagamma <- matrix(NA, nrow = breps, ncol = 6)

        for (b in 1:breps) {
          bcommunity <-
            community[sample(rownames(community), replace = TRUE),]
          bvar <- comm[rownames(bcommunity),sp]
          babun <- table(bvar)
          babun <- as.numeric(names(babun[babun >= 5]))

          for (k in babun[-1]) {
            b.deltagamma[b,k] <-   sum(colSums(bcommunity[bvar == k,]) > 0) - sum(colSums(bcommunity[bvar == 1,]) > 0)
          }
          rm(bcommunity, bvar)
          print(paste(
            sp, ": boot#", b, "gamma :",  paste(b.deltagamma[b,], collapse = " ")
          ))
        }
      }



      # NULL EXPECTATIONS
      r.gamma <- matrix(NA, nrow = nreps, ncol = 6)
      r.deltagamma <- matrix(NA, nrow = nreps, ncol = 6)

      for (j in 1:nreps) {
        rabun <- abun
        # null model A : simple permutation of all plots  ### keep sample size, not alphas, pool of species
        if (null.model == "permute.all") {
          rvar <- sample(var)
          rcommunity <- community
          rabun <- table(rvar)
          rabun <- as.numeric(names(rabun[rabun >= 5]))
        }

        # null model B1 : keep alpha richness, but resample sp. occurrences from total pool
        if( null.model == "resample.beta") {
          rcommunity <- community
          rvar <- var

          # for each site, resample identity of sp present, based on fequency of occurrence of sp across all plots
          rcommunity <-  t(apply(community, 1, FUN=function(x) {
            new <- rep(0,dim(community)[2] )
            new[ sample( 1:dim(community)[2],size = sum(x), prob = freq, replace =FALSE)] <- 1
            return(new)
          }))
          rfreq =colSums(rcommunity)/sum(rcommunity)
          rabun <- table(rvar)
          rabun <- as.numeric(names(rabun[rabun >= 5]))
        }

        # null model B2 : keep alpha richness, but resample sp. occurrences from total pool
        if( null.model == "resample.beta.inrare") {
          rcommunity <- community
          rvar <- var

          # for each site, resample identity of sp present, based on fequency of occurrence of sp across all plots
#       rcommunity <-  t(apply(community, 1, FUN=function(x) {
#             new <- rep(0,dim(community)[2] )
#             new[ sample( 1:dim(community)[2],size = sum(x), prob = freq.rare, replace =FALSE)] <- 1
#             return(new)
#           }))

          # for each site, resample identity of sp present, based on fequency of occurrence of sp across all plots
          rcommunity <-  t(sapply(rownames(community),FUN=function(x) {
            new <- rep(0,dim(community)[2] )
            new[ sample( 1:dim(community)[2],size = sum(community[x,]), prob = freq.rare, replace =FALSE)] <- 1
            return(new)
          }))
          rabun <- table(rvar)
          rabun <- as.numeric(names(rabun[rabun >= 5]))
        }

        ### calculate gamma diversity at each abundance
        for (k in rabun) {

          #Null model 0: permute with rare abundances
           if (null.model == "permute.rare") {
            rvar <- var
            rvar[rvar %in% c(1, k)] <-
              sample(rvar[rvar %in% c(1, k)])
            rcommunity <- community
          }

          # null model B : resample in Ab = (rare or K)
          if (null.model == "resample.beta.rarevsk") {
            rvar <- var
            # for each site, resample identity of sp present,
            # based on fequency of occurrence of sp across [abundance = rare or abundance = k plots]:
            subcom <- community[var %in% c(1, k),]
            subcom <- subcom[,colSums(subcom) > 0]

            rcommunity <-
              t(sapply(
                rowSums(community), FUN = function(x) {
                  new <- rep(0,dim(community)[2])
                  names(new) <- colnames(community)
                  new[sample(colnames(subcom),size = x, prob = colSums(subcom) / sum(var %in% c(1, k)), replace =
                      FALSE)] <- 1
                  return(new)
                }
              ))
          }

          r.gamma[j,k] <-  sum(colSums(rcommunity[rvar == k,]) > 0)
        }
        r.deltagamma[j,] <-r.gamma[j,] - r.gamma[j,1]
        print(r.gamma[j,])
        rm(rcommunity, rvar)
        print(paste(sp, ":" ,j))
      }
      r.gamma <-  rbind(gamma.obs$gamma,  r.gamma)
      r.deltagamma <-  rbind(gamma.obs$deltagamma, r.deltagamma)

      gamma.list$obs[sp,] <- gamma.obs$gamma
      gamma.list$alpha[sp,] <- gamma.obs$alpha
      gamma.list$P.gamma[sp,] <-
        apply(r.gamma, 2, function(x)
          (sum(x < x[1]) + sum(x == x[1]) / 2) / (nreps + 1))
      gamma.list$mean.null[sp,] <-
        apply(r.gamma, 2, function(x)
          mean(x, na.rm = T))
      gamma.list$sd.null[sp,] <-
        apply(r.gamma, 2, function(x)
          sd(x, na.rm = T))
      gamma.list$q975.null[sp,] <-
        apply(r.gamma, 2, function(x)
          quantile(x, probs = 0.975, na.rm = T))
      gamma.list$q025.null[sp,] <-
        apply(r.gamma, 2, function(x)
          quantile(x, probs = 0.025, na.rm = T))


      gamma.list$obs.deltagamma[sp,] <- gamma.obs$deltagamma
      gamma.list$P.deltagamma[sp,] <-
        apply(r.deltagamma, 2, function(x)
          (sum(x < x[1]) + sum(x == x[1]) / 2) / (nreps + 1))
      gamma.list$mean.null.deltagamma[sp,] <-
        apply(r.deltagamma, 2, function(x)
          mean(x, na.rm = T))
      gamma.list$sd.null.deltagamma[sp,] <-
        apply(r.deltagamma, 2, function(x)
          sd(x, na.rm = T))
      gamma.list$q975.null.deltagamma[sp,] <-
        apply(r.deltagamma, 2, function(x)
          quantile(x, probs = 0.975, na.rm = T))
      gamma.list$q025.null.deltagamma[sp,] <-
        apply(r.deltagamma, 2, function(x)
          quantile(x, probs = 0.025, na.rm = T))



      gamma.list$mean.boot[sp,] <-
        apply(b.deltagamma, 2, function(x)
          mean(x, na.rm = T))
      gamma.list$sd.boot[sp,] <-
        apply(b.deltagamma, 2, function(x)
          sd(x, na.rm = T))
      gamma.list$q975.boot[sp,] <-
        apply(b.deltagamma, 2, function(x)
          quantile(x, probs = 0.975, na.rm = T))
      gamma.list$q025.boot[sp,] <-
        apply(b.deltagamma, 2, function(x)
          quantile(x, probs = 0.025, na.rm = T))
    }

    return(gamma.list)
  }


## REwrite of div part funciton with all abundance threshold explorations
beta.trend <-
  function(spnames = impsp, group = natives, null.model = "permute.rare", nreps = 99) {
    # initiating output
    tmp <- as.data.frame(matrix(
      NA, nrow = length(spnames), ncol = 6,
      dimnames = list(spnames,paste("c", 1:6 , sep =
                                      ""))
    ))
    beta.list     <-
      list(
        obs = tmp , P.beta = tmp, mean.null = tmp, sd.null = tmp, q975 = tmp , q025 = tmp
      )
    rm(tmp)

    # looping on each of the target species
    for (sp in spnames) {
      community <-
        comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp] > 0)),]
      community <- community[,colSums(community) > 0]
      community <- community[,colnames(community) %in% group]
      community <- ceiling(community > 0)

      var <- comm[rownames(community),sp]

      freq <- colSums(community) / dim(community)[1]
      rare.freq <- colSums(community[var == 0,]) / sum(var == 0)

      divpart.obs <- c(
        n <- sum(rowSums(community[var == 1,]) > 0),
        gam <- sum(colSums(community[var == 1,]) > 0),
        alpha <-
          mean(rowSums(community[var == 1,]), na.rm = T),
        beta.add <- gam - alpha,
        beta.mult <- gam / alpha,
        beta.prop <- (gam - alpha) / gam
      )

      for (k in 2:6) {
        g <- (if (sum(var == k) > 1) {
          c(
            n <- sum(rowSums(community[var == k,]) > 0),
            gam <- sum(colSums(community[var == k,]) > 0),
            alpha <- mean(rowSums(community[var == k,]), na.rm = T),
            beta.add <- gam - alpha,
            beta.mult <- gam / alpha,
            beta.prop <- (gam - alpha) / gam
          )
        }
        else {
          c(
            n <- sum(sum(community[var == k,]) > 0),
            gam <- sum(community[var == k,] > 0),
            alpha <- sum(community[var == k,], na.rm = T),
            beta.add <- 0,
            beta.mult <- 0,
            beta.prop <- 0
          )
        })
        divpart.obs <- rbind(divpart.obs, g)
      }
      divpart.obs <- as.data.frame(divpart.obs)
      rownames(divpart.obs) <- paste("c", 1:6)
      colnames(divpart.obs) <-
        c("nb.plot", "gamma", "alpha", "beta.add", "beta.mult","beta.prop")


      # build random expectatinos
      r.betas <- matrix(NA, nrow = nreps, ncol = 6)

      for (j in 1:nreps) {
        # null model A : simple permutation of all plots
        if (null.model == "permute.all") {
          rvar <- sample(var)
        }

        ### calculate gamma diversity at each abundance
        r.betas[j,] <- sapply(
          1:6, FUN = function(k) {
            if (null.model == "permute.rare") {
              rvar <- var
              rvar[rvar %in% c(1, k)] <- sample(rvar[rvar %in% c(1, k)])
            }

            if (null.model == "permute.below") {
              rvar <- var
              rvar[rvar <= k] <- sample(rvar[rvar <= k])
            }

            if (sum(rvar == k) > 1) {
              # beta.prop
              b <-
                (sum(colSums(community[rvar == k,]) > 0) - mean(rowSums(community[var ==
                                                                                    k,]), na.rm = T)) / sum(colSums(community[rvar == k,]) > 0)

              #           #beta.miult
              #           b <- sum(colSums(community[rvar==k,])>0)/ mean(rowSums(community[var==k,]), na.rm=T)

            }
            if (sum(rvar == k) <= 1)
              b <- 0
            return(b)
          }
        )

        print(paste(sp, ":" ,j))
      }

      r.betas <-  rbind(divpart.obs$beta.prop,  r.betas)

      beta.list$obs[sp,] <- divpart.obs$beta.prop
      beta.list$P.beta[sp,] <-
        apply(r.betas, 2, function(x)
          (sum(x < x[1]) + sum(x == x[1]) / 2) / (nreps + 1))
      beta.list$mean.null[sp,] <-
        apply(r.betas, 2, function(x)
          mean(x, na.rm = T))
      beta.list$sd.null[sp,] <-
        apply(r.betas, 2, function(x)
          sd(x, na.rm = T))
      beta.list$q975[sp,] <-
        apply(r.betas, 2, function(x)
          quantile(x, probs = 0.975, na.rm = T))
      beta.list$q025[sp,] <-
        apply(r.betas, 2, function(x)
          quantile(x, probs = 0.025, na.rm = T))
    }

    return(beta.list)
  }
