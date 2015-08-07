### Partitioning Gamma diversity of communities
library(vegan)


## Overall partitionning of diversity
div.part.overall <-  function (){
community <- comm[which(rownames(comm) %in% realgrasslands),]
community <- community[,colSums(community)>0]
community<- ceiling(community>0)

n <- sum(rowSums(community)>0)
gam <- sum(colSums(community)>0)
alpha <- mean(rowSums(community), na.rm=T)
beta.add <- gam - alpha
beta.mult <- gam/alpha
beta.prop <- (gam - alpha)/gam

nbeta <- oecosimu(community, nestedbetasor, method = "r1", nsimul = 499, burnin = 0, thin = 1,
                        statistic = "statistic", alternative = c("less"))

turnover.obs <- nbeta$stat[1]
turnover.pval <- nbeta$oecosimu$pval[1]

nestedness.obs <- nbeta$stat[2]
nestedness.pval <- nbeta$oecosimu$pval[2]

div.part.overall <- c(n, gam, alpha, beta.add, beta.mult,beta.prop,
                      turnover.obs,turnover.pval,nestedness.obs,nestedness.pval)

comm.nat <- community[,colnames(community) %in% natives]

n <- sum(rowSums(comm.nat)>0)
gam <- sum(colSums(comm.nat)>0)
alpha <- mean(rowSums(comm.nat), na.rm=T)
beta.add <- gam - alpha
beta.mult <- gam/alpha
beta.prop <- (gam - alpha)/gam

nbeta <- oecosimu(community, nestedbetasor, method = "r1", nsimul = 499, burnin = 0, thin = 1,
                  statistic = "statistic", alternative = c("less"))

turnover.obs <- nbeta$stat[1]
turnover.pval <- nbeta$oecosimu$pval[1]

nestedness.obs <- nbeta$stat[2]
nestedness.pval <- nbeta$oecosimu$pval[2]

div.part.overall <- rbind(div.part.overall,
                          c(n, gam, alpha, beta.add, beta.mult,beta.prop,
                            turnover.obs,turnover.pval,nestedness.obs,nestedness.pval))

comm.ali <- community[,colnames(community) %in% aliens]
n <- sum(rowSums(comm.ali)>0)
gam <- sum(colSums(comm.ali)>0)
alpha <- mean(rowSums(comm.ali), na.rm=T)
beta.add <- gam - alpha
beta.mult <- gam/alpha
beta.prop <- (gam - alpha)/gam

nbeta <- oecosimu(community, nestedbetasor, method = "r1", nsimul = 499, burnin = 0, thin = 1,
                  statistic = "statistic", alternative = c("less"))

turnover.obs <- nbeta$stat[1]
turnover.pval <- nbeta$oecosimu$pval[1]

nestedness.obs <- nbeta$stat[2]
nestedness.pval <- nbeta$oecosimu$pval[2]

div.part.overall <- rbind(div.part.overall,
                          c(n, gam, alpha, beta.add, beta.mult,beta.prop,
                            turnover.obs,turnover.pval,nestedness.obs,nestedness.pval))

colnames(div.part.overall) = c("nb.plot", "gamm", "alpha", "beta.add", "beta.mult","beta.prop",
                              "turnover.obs","turnover.pval","nestedness.obs","nestedness.pval")
rownames(div.part.overall) = c("total richness", "native richness", "alien richness")
 return(div.part.overall)
}

## FUNCTION Partition diversity with standaridzation against null expectation
div.part.gamma.nm <- function(spnames = impsp, group = natives, thresh = c("th", "max"),
                        null.model = "permute", nreps = 99) {

  tmp<- as.data.frame(matrix(NA, nrow = length(spnames), ncol = 6,
                             dimnames= list(
                               spnames,
                              paste(rep("C"), 1:6, sep="")
                             )))

  divpart<- list( nplot.below = tmp,nplot.above = tmp,
                  alpha.below = tmp , alpha.above = tmp,
                  gamma.below = tmp , gamma.above = tmp,
                  betap.below = tmp , betap.above = tmp,
                  P.above = tmp, z.above= tmp, null.above = tmp, sdnull.above= tmp  )


  for (sp in spnames) {

    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
    community <- community[,colSums(community)>0]
    community<- community[,colnames(community) %in% group]
    community<- ceiling(community>0)

    var <- comm[rownames(community),sp]
    abun <- table(var)
    abun <- as.numeric(names(abun[abun>=5]))

    for (k in 1:max(abun)){
    divpart$alpha.below [sp,k] <- mean(rowSums(community[var<k,]), na.rm=T)
    divpart$alpha.above [sp,k] <- mean(rowSums(community[var>=k,]), na.rm=T)
    divpart$gamma.below [sp,k] <- sum(colSums(community[var<k,])>0)
    divpart$gamma.above [sp,k] <- sum(colSums(community[var>=k,])>0)
    divpart$nplot.below [sp,k] <- sum(var<k)
    divpart$nplot.above [sp,k] <- sum(var>=k)

    divpart$betap.below [sp,k] <- (divpart$gamma.below [sp,k] - divpart$alpha.below [sp,k]) /divpart$gamma.below [sp,k]
    divpart$betap.above [sp,k] <- (divpart$gamma.above [sp,k] - divpart$alpha.above [sp,k]) /divpart$gamma.above [sp,k]
  }


    # build random expectatinos
  r.gamma.above <- matrix(NA, nreps,6)

  for (j in 1:nreps) {

      # null model A : simple permutation of the below/above plots
      if( null.model == "permute") {
        rvar <- sample(var)
        rcommunity <- community
        abun <- table(rvar)
        abun <- as.numeric(names(abun[abun>=5]))
      }

      for (k in  1:max(abun)){
        r.gamma.above  [j,k] <- sum(colSums(rcommunity[rvar>=k,])>0)
             }
       print(paste(sp, ":" ,j))
    }
  r.gamma.above<-  rbind(as.numeric(divpart$gamma.above [sp,]), r.gamma.above)

    divpart$P.above[sp,] <- apply(r.gamma.above, 2, function(x) (sum(x<x[1]) + sum(x == x[1])/2)/(nreps+1))
    divpart$z.above[sp,]  <- apply(r.gamma.above, 2, function(x) {
      if(sd(x, na.rm=T) != 0 & !all(is.na(x))) z=(x[1] - mean(x, na.rm=T))/sd(x, na.rm=T)
      if(sd(x, na.rm=T) == 0 | all(is.na(x))) z=NA
      return(z)
    } )
  divpart$null.above[sp,]  <-  apply(r.gamma.above, 2, function(x) mean(x, na.rm=T))
  divpart$sdnull.above[sp,] <-  apply(r.gamma.above, 2, function(x) sd(x, na.rm=T))
  }

  return(divpart)
}

div.part.alpha.nm <- function(spnames = impsp, group = natives, thresh = c("th", "max"),
                        null.model = "permute", nreps = 99) {

  tmp<- as.data.frame(matrix(NA, nrow = length(spnames), ncol = 6,
                             dimnames= list(
                               spnames,
                               paste(rep("C"), 1:6, sep="")
                             )))

  divpart<- list( nplot.below = tmp,nplot.above = tmp,
                  alpha.below = tmp , alpha.above = tmp,
                  gamma.below = tmp , gamma.above = tmp,
                  betap.below = tmp , betap.above = tmp,
                  P.above = tmp, z.above= tmp, null.above = tmp, sdnull.above= tmp  )


  for (sp in spnames) {

    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
    community <- community[,colSums(community)>0]
    community<- community[,colnames(community) %in% group]
    community<- ceiling(community>0)

    var <- comm[rownames(community),sp]
    abun <- table(var)
    abun <- as.numeric(names(abun[abun>=5]))

    for (k in 1:max(abun)){
      divpart$alpha.below [sp,k] <- mean(rowSums(community[var<k,]), na.rm=T)
      divpart$alpha.above [sp,k] <- mean(rowSums(community[var>=k,]), na.rm=T)
      divpart$gamma.below [sp,k] <- sum(colSums(community[var<k,])>0)
      divpart$gamma.above [sp,k] <- sum(colSums(community[var>=k,])>0)
      divpart$nplot.below [sp,k] <- sum(var<k)
      divpart$nplot.above [sp,k] <- sum(var>=k)

      divpart$betap.below [sp,k] <- (divpart$gamma.below [sp,k] - divpart$alpha.below [sp,k]) /divpart$gamma.below [sp,k]
      divpart$betap.above [sp,k] <- (divpart$gamma.above [sp,k] - divpart$alpha.above [sp,k]) /divpart$gamma.above [sp,k]
    }


    # build random expectatinos
    r.alpha.above <- matrix(NA, nreps,6)

    for (j in 1:nreps) {

      # null model A : simple permutation of the below/above plots
      if( null.model == "permute") {
        rvar <- sample(var)
        rcommunity <- community
        abun <- table(rvar)
        abun <- as.numeric(names(abun[abun>=5]))
      }

      for (k in  1:max(abun)){
        r.alpha.above  [j,k] <- mean(rowSums(community[rvar>=k,]), na.rm=T)
      }
      print(paste(sp, ":" ,j))
    }
    r.alpha.above<-  rbind(as.numeric(divpart$alpha.above [sp,]), r.alpha.above)

    divpart$P.above[sp,] <- apply(r.alpha.above, 2, function(x) (sum(x<x[1]) + sum(x == x[1])/2)/(nreps+1))
    divpart$z.above[sp,]  <- apply(r.alpha.above, 2, function(x) {
      if(sd(x, na.rm=T) != 0 & !all(is.na(x))) z=(x[1] - mean(x, na.rm=T))/sd(x, na.rm=T)
      if(sd(x, na.rm=T) == 0 | all(is.na(x))) z=NA
      return(z)
    } )
    divpart$null.above[sp,]  <-  apply(r.alpha.above, 2, function(x) mean(x, na.rm=T))
    divpart$sdnull.above[sp,] <-  apply(r.alpha.above, 2, function(x) sd(x, na.rm=T))
  }

  return(divpart)
}




