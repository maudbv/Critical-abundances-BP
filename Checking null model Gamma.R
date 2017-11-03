# CHecking that the null model for gamma richness is not too biased by small sample sizes.

## number of potential null expectations for p plots sampled from at least n other plots
# n!/(n! * (n-p)!) = Combinaison de p élements distincts parmis n

# number of unique combinations of 5 plots from a pool of ten plots:
factorial(10)/(factorial(5) * factorial(5)) # only 252

# number of unique combinations of 5 plots from a pool of ten plots: (minimal scenario)
factorial(10)/(factorial(5) * factorial(5)) # only 252

# number of unique combinations of 6 plots from a pool of 52 plots:(Achillea millefolium)
factorial(52)/(factorial(6) * factorial(52-6)) # 2598960

# number of unique combinations of 6 plots from a pool of 82 plots:(Anthoxanthum odoratum)
factorial(82)/(factorial(5) * factorial(82-5)) # 27285336


# Present histograms of null expectations
spnames = impsp
group = natives
bootstrapped = TRUE
null.model = "permute.rare"
breps= 99
nreps = 999

  # looping on each of the target species
sp <-  spnames[1]
sp <- "POACIT"
    community <-comm[which((rownames(comm) %in% unimprovedgrasslands) & (comm[,sp] > 0)),]
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
        n <- sum(var == k),
        gam <- sum(colSums(community[var == k,]) > 0),
        alpha <- mean(rowSums(community[var == k,]), na.rm = T),
        beta.add <- gam - alpha,
        beta.mult <- gam / alpha,
        beta.prop <- (gam - alpha) / gam
      )
    }
    gamma.obs $deltagamma <-gamma.obs$gamma - gamma.obs$gamma[1]
    
    # BOOTSTRAPPING delta gamma
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
    
      
      par(mfrow = c(1,4))
      tmp <- rbind(gamma.obs[, "deltagamma"], b.deltagamma)
      hist(tmp[,2], col = "grey")
      abline(v = tmp[1,2], col = "red")
      hist(tmp[,3], col = "grey")
      abline(v =tmp[1,3], col = "red")
      hist(tmp[,4], col = "grey")
      abline(v = tmp[1,4], col = "red")
      hist(tmp[,5], col = "grey")
      abline(v = tmp[1,5], col = "red")
      
      par(mfrow = c(1,1))
      plot(tmp[1,], col = "red", type = "b", ylim = c(min(tmp, na.rm = T), max(tmp, na.rm  =T)))
      for (i in 1:breps+1) lines(tmp[i,], type = "b", col = "grey")
      lines(tmp[1,], col = "red", type = "b")
    
    # NULL EXPECTATIONS
    r.gamma <- matrix(NA, nrow = nreps, ncol = 6)
    r.deltagamma <- matrix(NA, nrow = nreps, ncol = 6)
    
    for (j in 1:nreps) {
      rabun <- abun
      ### calculate gamma diversity at each abundance
      for (k in rabun) {
        #Null model 0: permute with rare abundances
        if (null.model == "permute.rare") {
          rvar <- var
          rvar[rvar %in% c(1, k)] <-
            sample(rvar[rvar %in% c(1, k)])
          rcommunity <- community
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
    
    
    par(mfrow = c(1,5))
    tmp <- rbind(gamma.obs[, "gamma"], r.gamma)
    hist(tmp[,2], col = "grey", main = "class2", xlab = "Native Gamma richness")
    abline(v = tmp[1,2], col = "red")
    hist(tmp[,3], col = "grey", main = "class3", xlab = "Native Gamma richness")
    abline(v =tmp[1,3], col = "red")
    hist(tmp[,4], col = "grey", main = "class4", xlab = "Native Gamma richness")
    abline(v = tmp[1,4], col = "red")
    hist(tmp[,5], col = "grey", main = "class5", xlab = "Native Gamma richness")
    abline(v = tmp[1,5], col = "red")
    hist(tmp[,6], col = "grey", main = "class6", xlab = "Native Gamma richness")
    abline(v = tmp[1,6], col = "red")
    
    par(mfrow = c(1,1))
    plot(tmp[1,], col = "red", type = "n", ylim = c(min(tmp, na.rm = T), max(tmp, na.rm  =T)),
         main = sp, ylab = "Native Gamma richness", xlab = "abundance class" )
    for (i in 1:breps+1) lines(tmp[i,], type = "b", col = "grey")
    lines(tmp[1,], col = "red", type = "b")
    mtext(3, line = 0, text  =  gamma.obs[, "nb.plot"], at = 1:6, col = "darkgrey")
    
