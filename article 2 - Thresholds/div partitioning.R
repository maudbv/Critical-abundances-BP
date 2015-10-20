### Partitioning Gamma diversity of communities above and below thresholds
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

#### randomize alphas : shuffle sites above/below
div.part.alpha.nm <- function(spnames = impsp, group = natives, thresh = c("th", "max"),
                        null.model = "permute", nreps = 99) {

  tmp<- as.data.frame(matrix(NA, nrow = length(spnames), ncol = 6,
                             dimnames= list(
                               spnames,
                               paste(rep("C"), 1:6, sep="")
                             )))

  divpart<- list( nplot.below = tmp,nplot.above = tmp,
                  alpha.below = tmp , alpha.above = tmp,
                  alpha.below.sd = tmp , alpha.above.sd = tmp,
                  deltaalpha = tmp ,
                  gamma.below = tmp , gamma.above = tmp,
                  betap.below = tmp , betap.above = tmp,
                  t=tmp, df = tmp, P.t = tmp,
                  P.above = tmp, z.above= tmp, null.above = tmp, sdnull.above= tmp,
                  null.delta = tmp, sdnull.delta= tmp)


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
      divpart$alpha.below.sd [sp,k] <- sd(rowSums(community[var<k,]), na.rm=T)
      divpart$alpha.above.sd [sp,k] <- sd(rowSums(community[var>=k,]), na.rm=T)
      divpart$deltaalpha [sp,k] <-divpart$alpha.above.sd [sp,k] - divpart$alpha.below.sd [sp,k]
      divpart$gamma.below [sp,k] <- sum(colSums(community[var<k,])>0)
      divpart$gamma.above [sp,k] <- sum(colSums(community[var>=k,])>0)
       divpart$nplot.below [sp,k] <- sum(var<k)
      divpart$nplot.above [sp,k] <- sum(var>=k)

      t.above = list( NA,NA,NA)
      if (sum(var<k)>1 & sum(var>=k)>1 ) t.above <- t.test(rowSums(community[var<k,]), rowSums(community[var>=k,]))

      divpart$t [sp,k] <- t.above[1]
      divpart$df[sp,k] <- t.above[2]
      divpart$P.t[sp,k] <- t.above[3]

      divpart$betap.below [sp,k] <- (divpart$gamma.below [sp,k] - divpart$alpha.below [sp,k]) /divpart$gamma.below [sp,k]
      divpart$betap.above [sp,k] <- (divpart$gamma.above [sp,k] - divpart$alpha.above [sp,k]) /divpart$gamma.above [sp,k]
    }


    # build random expectatinos
    r.alpha.above <- matrix(NA, nreps,6)
    r.deltaalpha <- matrix(NA, nreps,6)
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
        r.deltaalpha  [j,k] <- mean(rowSums(community[rvar>=k,]), na.rm=T) - mean(rowSums(community[rvar<k,]), na.rm=T)
      }
      print(paste(sp, ":" ,j))
    }
    r.alpha.above<-  rbind(as.numeric(divpart$alpha.above [sp,]), r.alpha.above)
    r.deltaalpha<-  rbind(as.numeric(divpart$deltaapha [sp,]), r.deltaalpha)

    divpart$P.above[sp,] <- apply(r.alpha.above, 2, function(x) (sum(x<x[1]) + sum(x == x[1])/2)/(nreps+1))
    divpart$z.above[sp,]  <- apply(r.alpha.above, 2, function(x) {
      if(sd(x, na.rm=T) != 0 & !all(is.na(x))) z=(x[1] - mean(x, na.rm=T))/sd(x, na.rm=T)
      if(sd(x, na.rm=T) == 0 | all(is.na(x))) z=NA
      return(z)
    } )
    divpart$null.above[sp,]  <-  apply(r.alpha.above, 2, function(x) mean(x, na.rm=T))
    divpart$sdnull.above[sp,] <-  apply(r.alpha.above, 2, function(x) sd(x, na.rm=T))

    divpart$null.delta[sp,] <- apply( r.deltaalpha, 2, function(x) mean(x, na.rm=T))
    divpart$sdnull.delta[sp,] <- apply( r.deltaalpha, 2, function(x) sd(x, na.rm=T))


  }

  return(divpart)
}



#  using betadisper test by Anderson et al. to look for drift in betadiersity above threshold
# betasim.test <- data.frame(matrix(NA, nrow = length(spnames), ncol = 6, dimnames = list(spnames, 1:6)))
# par (mfrow = c(3,4), mar = c(2,2,2,1))
#   for (sp in spnames) {
#
#     community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
#     community <- community[,colSums(community)>0]
#     community<- community[,colnames(community) %in% group]
#     community<- ceiling(community>0)
#     community <- community[rowSums(community) > 0 , ] ### have to remove empty communities
#
#     var <- comm[rownames(community),sp]
#     abun <- table(var)
#     abun <- as.numeric(names(abun[abun>=5]))
#
#     g<- var %in% abun
#
#     for (k in abun [-1]) {
#     vark <- as.numeric(var[g]>=k)
#     sim.dist <- betadiver(community[g,], method = "sim")   # pmin(b+c) / (a + pmin(b + c))
#     mod <- betadisper(sim.dist, group = vark)
#     p.test <- permutest(mod, permutations = 999)
#     betasim.test[sp,k] <- p.test$tab$`Pr(>F)`[1]
#   }
# }



## FUNCTION Model loss in gamma diversity, bootstrapped and with standaridzation against null expectation
div.part.gamma.nm <- function(spnames = impsp, group = natives, thresh = c("th", "max"),
                              null.model = "permute", nreps = 99,
                              bootstrapped = TRUE, breps = 99) {

  tmp<- as.data.frame(matrix(NA, nrow = length(spnames), ncol = 6,
                             dimnames= list(
                               spnames,
                               paste(rep("C"), 1:6, sep="")
                             )))

  divpart<- list( nplot.below = tmp,nplot.above = tmp,
                  alpha.below = tmp , alpha.above = tmp,
                  gamma.below = tmp , gamma.above = tmp,
                  deltagamma = tmp ,
                  betap.below = tmp , betap.above = tmp,
                  P.above = tmp, z.above= tmp, null.above = tmp, sdnull.above= tmp,
                  P.delta = tmp, z.delta= tmp, null.delta = tmp, sdnull.delta= tmp,
                  lowCI.null.above = tmp, hiCI.null.above = tmp,
                  lowCI.null.delta = tmp, hiCI.null.delta = tmp,
                  lowCI.boot.delta = tmp, hiCI.boot.delta= tmp, mean.boot.delta= tmp
                  )


  for (sp in spnames) {

    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
    community <- community[,colSums(community)>0]
    community<- community[,colnames(community) %in% group]
    community<- ceiling(community>0)

    var <- comm[rownames(community),sp]
    abun <- table(var)
    abun <- as.numeric(names(abun[abun>=5]))

    freq <- colSums(community) / dim(community)[1]
    freq.rare <- colSums(community[var==1,]) / dim(community[var==1,])[1]

    for (k in 1:max(abun)){
      divpart$alpha.below [sp,k] <- mean(rowSums(community[var<k,]), na.rm=T)
      divpart$alpha.above [sp,k] <- mean(rowSums(community[var>=k,]), na.rm=T)
      divpart$gamma.below [sp,k] <- sum(colSums(community[var<k,])>0)
      divpart$gamma.above [sp,k] <- sum(colSums(community[var>=k,])>0)
      divpart$deltagamma [sp,k] <- sum(colSums(community[var>=k,])>0) - sum(colSums(community[var<k,])>0)

      divpart$nplot.below [sp,k] <- sum(var<k)
      divpart$nplot.above [sp,k] <- sum(var>=k)

      divpart$betap.below [sp,k] <- (divpart$gamma.below [sp,k] - divpart$alpha.below [sp,k]) /divpart$gamma.below [sp,k]
      divpart$betap.above [sp,k] <- (divpart$gamma.above [sp,k] - divpart$alpha.above [sp,k]) /divpart$gamma.above [sp,k]
    }


    # BOOTSTRAPPING gamma above values
    if (bootstrapped) {
      # build random expectations
      b.deltagamma <- matrix(NA, nrow = breps, ncol = 6)
      for (b in 1:breps) {
        bcommunity <-
        community[sample(rownames(community), replace = TRUE),]
        bvar <- comm[rownames(bcommunity),sp]
        babun <- table(bvar)
        babun <- as.numeric(names(babun[babun >= 5]))

        for (k in babun) {
        b.deltagamma[b,k] <-   sum(colSums(bcommunity[bvar>=k,])>0) - sum(colSums(bcommunity[bvar<k,])>0)
        }
        rm(bcommunity, bvar)
        print(paste(
          sp, ": boot#", b, "gamma :",  paste(b.deltagamma[b,], collapse = " ")
        ))
      }
        b.deltagamma<-  rbind(as.numeric(divpart$deltagamma [sp,]), b.deltagamma)

     # store bootstrap results
     divpart$lowCI.boot.delta[sp,]  <-  apply( b.deltagamma, 2, function(x) quantile(x, 0.025,na.rm=T))
     divpart$hiCI.boot.delta[sp,] <-  apply( b.deltagamma, 2, function(x)  quantile(x, 0.975,na.rm=T))
     divpart$mean.boot.delta[sp,]  <-  apply( b.deltagamma, 2, function(x) mean(x, na.rm=T))
    }

    # NULL MODELLING
    r.gamma.above <- matrix(NA, nreps,6)
    r.deltagamma <- matrix(NA, nreps,6)

    for (j in 1:nreps) {

      # null model 0a : simple permutation of the plots across all grassland plots
      if( null.model == "permute.total") {
        rcommunity <- comm[(rownames(comm) %in% realgrasslands),]
        rcommunity <- rcommunity[,colSums(rcommunity)>0]
        rcommunity<- rcommunity[,colnames(rcommunity) %in% group]
        rcommunity<- ceiling(rcommunity>0)

        rvar <- var
        rabun <- abun
        rcommunity <- rcommunity[sample(rownames(rcommunity),length(var)),]

      }
      # null model 0 : simple permutation of the plots
      if( null.model == "permute.all") {
        rvar <- sample(var)
        rcommunity <- community
        rabun <- abun
      }

      # null model 1 : keep alpha richness distribution above/below, but change species distrib (= beta div)
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

      # NM2 : Resample in rare
      if( null.model == "resample.beta.inrare") {
        rvar <- var
        # for each site, resample identity of sp present, based on fequency of occurrence of sp across all plots
        rcommunity <-  t(sapply(rownames(community),FUN=function(x) {
          new <- rep(0,dim(community)[2] )
          new[ sample( 1:dim(community)[2],size = sum(community[x,]), prob = freq.rare, replace =FALSE)] <- 1
          return(new)
        }))
        rabun <- table(rvar)
        rabun <- as.numeric(names(rabun[rabun >= 5]))
      }


      for (k in  1:max(rabun)){
        r.gamma.above  [j,k] <- sum(colSums(rcommunity[rvar>=k,])>0)
        r.deltagamma[j,k] <-   sum(colSums(rcommunity[rvar>=k,])>0) - sum(colSums(rcommunity[rvar<k,])>0)
      }
      print(paste(sp, ":" ,j))
      remove(rcommunity, rvar)
    }
    r.gamma.above<-  rbind(as.numeric(divpart$gamma.above [sp,]), r.gamma.above)
    r.deltagamma<-  rbind(as.numeric(divpart$deltagamma [sp,]), r.deltagamma)


    # NM results

    # Test gamma above
    divpart$P.above[sp,] <- apply(r.gamma.above, 2, function(x) (sum(x<x[1]) + sum(x == x[1])/2)/(nreps+1))
    divpart$z.above[sp,]  <- apply(r.gamma.above, 2, function(x) {
      z <-(  if(!is.na(sd(x, na.rm=T)) & sd(x, na.rm=T) != 0 & !all(is.na(x))) z=(x[1] - mean(x, na.rm=T))/sd(x, na.rm=T)
             else z=NA )
         return(z)
    } )

    divpart$null.above[sp,]  <-  apply(r.gamma.above, 2, function(x) mean(x, na.rm=T))
    divpart$sdnull.above[sp,] <-  apply(r.gamma.above, 2, function(x) sd(x, na.rm=T))
    divpart$lowCI.null.above[sp,]  <-  apply( r.gamma.above, 2, function(x) quantile(x, 0.025,na.rm=T))
    divpart$hiCI.null.above[sp,] <-  apply( r.gamma.above, 2, function(x)  quantile(x, 0.975,na.rm=T))


    #Test Delta Gamma
    divpart$P.delta[sp,] <- apply(r.deltagamma, 2, function(x) (sum(x<x[1]) + sum(x == x[1])/2)/(nreps+1))
    divpart$z.delta[sp,]  <- apply(r.deltagamma, 2, function(x) {
      z <-(  if(!is.na(sd(x, na.rm=T)) & sd(x, na.rm=T) != 0 & !all(is.na(x))) z=(x[1] - mean(x, na.rm=T))/sd(x, na.rm=T)
             else z=NA )
      return(z)
    } )
    divpart$null.delta[sp,]  <-  apply(r.deltagamma, 2, function(x) mean(x, na.rm=T))
    divpart$sdnull.delta[sp,] <-  apply(r.deltagamma, 2, function(x) sd(x, na.rm=T))
    divpart$lowCI.null.delta[sp,]  <-  apply( r.deltagamma, 2, function(x) quantile(x, 0.025,na.rm=T))
    divpart$hiCI.null.delta[sp,] <-  apply( r.deltagamma, 2, function(x)  quantile(x, 0.975,na.rm=T))


    # Boostrap CI for delta gamma
    divpart$lowCI.boot.delta[sp,]  <-  apply( b.deltagamma, 2, function(x) quantile(x, 0.025,na.rm=T))
    divpart$hiCI.boot.delta[sp,] <-  apply( b.deltagamma, 2, function(x)  quantile(x, 0.975,na.rm=T))
    divpart$mean.boot.delta[sp,]  <-  apply(b.deltagamma, 2, function(x) mean(x, na.rm=T))

    }

  return(divpart)
}

