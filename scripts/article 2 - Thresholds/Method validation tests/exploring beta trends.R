### Variations in BETA diversity at each abundance for all species 
# (or at least the "important species" from impsp to start with)

## calculate beta (Simpson) for each abundance of a species:
beta.prop <-  function(sp = "ACHMIL",  com = community, group = natives)  {
  com <- com[which((rownames(com) %in% realgrasslands)), ]
  com <- com[, colSums(com) > 0]
  var <- com[, which(colnames(com) == sp)]
  com <- com[, colnames(com) %in% group]
  com <- ceiling(com > 0)

  betas <- sapply(1:6, function(k) {
    if (sum(var == k) > 1) {
      beta = beta.multi(com[var == k, ])$beta.SIM
    }
    else {
      beta = NA
    }
    return(beta)
  })
    return(betas)
}

#for all focal species:
comm.betas.native <- sapply(rownames(glmSRnat.overall$glms), beta.prop, group = natives)
comm.betas.aliens <- sapply(rownames(glmSRnat.overall$glms), beta.prop, group = aliens)

# compare the trends:
plot(1:6, rep(0:1, 3), type = "n", ann = F)
apply(comm.betas.natives , 2 , function(y) points(1:6, y))
apply(comm.betas.natives[, impsp] , 2 , function(y) points(1:6, y, col = "red", type = "b"))

# compare the trends:
plot(1:6, rep(0:1, 3), type = "n", ann = F)
apply(comm.betas.aliens , 2 , function(y) points(1:6, y))
apply(comm.betas.aliens[, impsp] , 2 , function(y) points(1:6, y, col = "red", type = "b"))

### bootstrapping for a species
boot.betas <- sapply(1:99, function(rep) {
  b <- beta.prop(sp = "ACHMIL", com = community[sample(nrow(community), replace = T),])
  return(b)
})



# build random expectatinos

r.betas <- matrix(NA, nrow = nreps, ncol = 6)

for (j in 1:nreps) {
  
  # null model A : simple permutation of all plots
  if (null.model == "permute.all") {
    
  rbetas <- t(sapply(spnames, function(sp)  {
    community <- comm[which((rownames(comm) %in% realgrasslands)), ]
    community <- community[, colSums(community) > 0]
    community <- community[, colnames(community) %in% group]
    community <- ceiling(community > 0)
    
    var <- comm[rownames(community), sp]
    freq <- colSums(community) / dim(community)[1]
    rare.freq <- colSums(community[var == 0, ]) / sum(var == 0)
    
    # random samples of plots
    rvar <- sample(var)
    
    # calculate beta metrics
    betas <- sapply(1:6, function(k) {
      if (sum(var == k) > 1) {
        n <- sum(rowSums(community[var == k, ]) > 0)
        gam <- sum(colSums(community[var == k, ]) > 0)
        alpha <- mean(rowSums(community[var == k, ]), na.rm = T)
        beta  <-  (gam - alpha) / gam  ## proportional beta
      }
      else {
        beta = 0
      }
      return(beta)
    })
    return(betas)
  }, simplify = T))

}

print(paste(sp, ":" , j))

r.betas <-  rbind(divpart.obs$beta.prop,  r.betas)

beta.list$obs[sp, ] <- divpart.obs$beta.prop
beta.list$P.beta[sp, ] <-
  apply(r.betas, 2, function(x)
    (sum(x < x[1]) + sum(x == x[1]) / 2) / (nreps + 1))
beta.list$mean.null[sp, ] <-
  apply(r.betas, 2, function(x)
    mean(x, na.rm = T))
beta.list$sd.null[sp, ] <-
  apply(r.betas, 2, function(x)
    sd(x, na.rm = T))
beta.list$q975[sp, ] <-
  apply(r.betas, 2, function(x)
    quantile(x, probs = 0.975, na.rm = T))
beta.list$q025[sp, ] <-
  apply(r.betas, 2, function(x)
    quantile(x, probs = 0.025, na.rm = T))
}




## REwrite of div part funciton with all abundance threshold explorations
beta.trend <-
  function(spnames = impsp,
           group = natives,
           null.model = "permute.rare",
           nreps = 99) {
    # initiating output
    tmp <- as.data.frame(matrix(
      NA,
      nrow = length(spnames),
      ncol = 6,
      dimnames = list(spnames, paste("c", 1:6 , sep =
                                       ""))
    ))
    beta.list     <-
      list(
        obs = tmp ,
        P.beta = tmp,
        mean.null = tmp,
        sd.null = tmp,
        q975 = tmp ,
        q025 = tmp
      )
    rm(tmp)
    
    # looping on each of the target species
    for (sp in spnames) {
      community <-
        comm[which((rownames(comm) %in% realgrasslands) &
                     (comm[, sp] > 0)), ]
      community <- community[, colSums(community) > 0]
      community <- community[, colnames(community) %in% group]
      community <- ceiling(community > 0)
      
      var <- comm[rownames(community), sp]
      
      freq <- colSums(community) / dim(community)[1]
      rare.freq <- colSums(community[var == 0, ]) / sum(var == 0)
      
      divpart.obs <- c(
        n <- sum(rowSums(community[var == 1, ]) > 0),
        gam <- sum(colSums(community[var == 1, ]) > 0),
        alpha <-
          mean(rowSums(community[var == 1, ]), na.rm = T),
        beta.add <- gam - alpha,
        beta.mult <- gam / alpha,
        beta.prop <- (gam - alpha) / gam
      )
      
      for (k in 2:6) {
        g <- (if (sum(var == k) > 1) {
          c(
            n <- sum(rowSums(community[var == k, ]) > 0),
            gam <- sum(colSums(community[var == k, ]) > 0),
            alpha <- mean(rowSums(community[var == k, ]), na.rm = T),
            beta.add <- gam - alpha,
            beta.mult <- gam / alpha,
            beta.prop <- (gam - alpha) / gam
          )
        }
        else {
          c(
            n <- sum(sum(community[var == k, ]) > 0),
            gam <- sum(community[var == k, ] > 0),
            alpha <- sum(community[var == k, ], na.rm = T),
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
        c("nb.plot",
          "gamma",
          "alpha",
          "beta.add",
          "beta.mult",
          "beta.prop")
      
      
      # build random expectatinos
      r.betas <- matrix(NA, nrow = nreps, ncol = 6)
      
      for (j in 1:nreps) {
        # null model A : simple permutation of all plots
        if (null.model == "permute.all") {
          rvar <- sample(var)
        }
        
        ### calculate gamma diversity at each abundance
        r.betas[j, ] <- sapply(
          1:6,
          FUN = function(k) {
            if (null.model == "permute.rare") {
              rvar <- var
              rvar[rvar %in% c(1, k)] <-
                sample(rvar[rvar %in% c(1, k)])
            }
            
            if (null.model == "permute.below") {
              rvar <- var
              rvar[rvar <= k] <- sample(rvar[rvar <= k])
            }
            
            if (sum(rvar == k) > 1) {
              # beta.prop
              b <-
                (sum(colSums(community[rvar == k, ]) > 0) - mean(rowSums(community[var ==
                                                                                     k, ]), na.rm = T)) / sum(colSums(community[rvar == k, ]) > 0)
              
              #           #beta.miult
              #           b <- sum(colSums(community[rvar==k,])>0)/ mean(rowSums(community[var==k,]), na.rm=T)
              
            }
            if (sum(rvar == k) <= 1)
              b <- 0
            return(b)
          }
        )
        
        print(paste(sp, ":" , j))
      }
      
      r.betas <-  rbind(divpart.obs$beta.prop,  r.betas)
      
      beta.list$obs[sp, ] <- divpart.obs$beta.prop
      beta.list$P.beta[sp, ] <-
        apply(r.betas, 2, function(x)
          (sum(x < x[1]) + sum(x == x[1]) / 2) / (nreps + 1))
      beta.list$mean.null[sp, ] <-
        apply(r.betas, 2, function(x)
          mean(x, na.rm = T))
      beta.list$sd.null[sp, ] <-
        apply(r.betas, 2, function(x)
          sd(x, na.rm = T))
      beta.list$q975[sp, ] <-
        apply(r.betas, 2, function(x)
          quantile(x, probs = 0.975, na.rm = T))
      beta.list$q025[sp, ] <-
        apply(r.betas, 2, function(x)
          quantile(x, probs = 0.025, na.rm = T))
    }
    
    return(beta.list)
  }
