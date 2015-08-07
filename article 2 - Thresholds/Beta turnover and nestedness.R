
# Beta diversity partitioning
betasor<- function(spnames = impsp, group = natives, nreps = 99) {

  betasor.list    <- list()

  for (i in 1:length(spnames)) {
    sp = spnames[i]
    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
    community <- community[,colSums(community)>0]
    community<- community[,colnames(community) %in% group]
    community<- ceiling(community>0)

    var <- comm[rownames(community),sp]
    var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
    var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1


    # beta partitioning
    nbeta.below <- oecosimu(community[var==0,], nestedbetasor, method = "r1", nsimul = nreps, burnin = 0, thin = 1,
                            statistic = "statistic", alternative = c("less"))

    nbeta.above <- oecosimu(community[var==1,], nestedbetasor, method = "r1", nsimul = nreps, burnin = 0, thin = 1,
                            statistic = "statistic", alternative ="less")

    # nestbetasor.SRali <- oecosimu(comm.ali, nestedbetasor, method = "r1", nsimul = 99, burnin = 0, thin = 1,
    #                               statistic = "statistic", alternative = c("two.sided", "less", "greater"))

    betasor.list[[i]] <- cbind(below = data.frame(obs = nbeta.below$oecosimu$statistic,
                                                  mean=nbeta.below$oecosimu$means,
                                                  z=nbeta.below$oecosimu$z,
                                                  pval = nbeta.below$oecosimu$pval),
                               above = data.frame(obs = nbeta.above$oecosimu$statistic,
                                                  mean=nbeta.above$oecosimu$means,
                                                  z=nbeta.above$oecosimu$z,
                                                  pval = nbeta.above$oecosimu$pval))


    names(betasor.list)[i] <- sp

  }
  return(betasor.list)

}
