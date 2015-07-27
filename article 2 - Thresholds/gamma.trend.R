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






###
# #########    plot GAMMA diversity  ##########
# 
# x11()
# par(mfrow = c(4,3), mar=c(4,4,1,1))
# for (i in 1 : length(impsp)){
# 
# plot(c(1,6.3),range(gamma.trend.nat$obs[i,],gamma.trend.nat.permute$q975[i,]), type = "n", ann=F, las = 1)
# segments(1:6,t(gamma.trend.nat$q025[i,]),
#          1:6,t(gamma.trend.nat$q975[i,]), col="grey", lwd = 3)
# 
# segments(1.1:6.1,t(gamma.trend.nat.below$q025[i,]),
#          1.1:6.1,t(gamma.trend.nat.below$q975[i,]), col="palegreen", lwd = 3)
# segments(1.2:6.2,t(gamma.trend.nat.below$q025[i,]),
#          1.2:6.2,t(gamma.trend.nat.below$q975[i,]), col="goldenrod", lwd = 3)
# 
# points(1:6,gamma.trend.nat$mean[i,],pch= 20, col="darkgrey")
# points(1:6,gamma.trend.nat$obs[i,])
# mtext(3, text = sub("_", " ",species[impsp[i],"tip"]), cex = 0.8)
# }
# 
# plot(1:10,1:10,type = "n", ann=F,axes=F)
# legend('left', legend = c("observed", "permute with rare occurrences", 
#                             "permute with lower abundances","permute with all plots"),
#        lty=c(NA, rep("solid",3)), lwd=c(0,5,5,5), pch=c(21, NA, NA, NA),
#              col=c("black", "grey", "palegreen", "goldenrod"), bty= "n",
#        cex = 0.8)
# 
# ## plot Standardized effect sizes
# x11()
# par(mfrow = c(4,3), mar=c(4,4,1,1))
# for (i in 1 : length(impsp)){
# cols <- c(NA, "black") [ (gamma.trend.nat$P.gamma[i,] <=0.025) +1]
# x = 1:6
# y =(gamma.trend.nat$obs[i,] - gamma.trend.nat$mean[i,])/ gamma.trend.nat$sd[i,]
# y[1]=0
# plot(x,y,ylim= c(-5,4),xlim= c(1,6.3),ann=F, las = 1, pch = 21, col ="black",bg = cols, type= "b")
# 
# 
# cols <- c(NA, "forestgreen") [ (gamma.trend.nat.below$P.gamma[i,] <=0.025) +1]
# x = 1.1:6.1
# y =(gamma.trend.nat.below$obs[i,] - gamma.trend.nat.below$mean[i,])/ gamma.trend.nat.below$sd[i,]
# par(new=T)
# plot(x,y,ylim= c(-5,4),xlim= c(1,6.3),pch = 21, col ="forestgreen",bg = cols, type= "b",ann=F, axes=F)
# 
# cols <- c(NA, "goldenrod") [ (gamma.trend.nat.permute$P.gamma[i,] <=0.025) +1]
# x = 1.2:6.2
# y =(gamma.trend.nat.permute$obs[i,] - gamma.trend.nat.permute$mean[i,])/ gamma.trend.nat.permute$sd[i,]
# par(new=T)
# plot(x,y,ylim= c(-5,4),xlim= c(1,6.3),pch = 21, col ="goldenrod",bg = cols, type= "b",ann=F, axes = F)
# 
# abline(h=0, lty="dotted", col="darkgrey")
# mtext(3, text = impsp[i], cex = 0.8)
# }
# plot(1:10,1:10,type = "n", ann=F,axes=F)
# legend('left', legend = c("permute with rare occurrences", 
#                           "permute with lower abundances","permute with all plots"),
#        lty=rep("solid",3),  pch=c(21, 21, 21),
#        col=c("black", "palegreen", "goldenrod"),
#        bty= "n",cex = 0.8)
# 
# 
# ########### plot BETA diversity #############
# 
# 
# x11()
# par(mfrow = c(4,3), mar=c(4,4,1,1))
# for (i in 1 : length(impsp)){
#   
#   plot(c(1,6.3),range(beta.trend.nat$obs[i,],beta.trend.nat.permute$q975[i,], na.rm = T), type = "n", ann=F, las = 1)
#   segments(1:6,t(beta.trend.nat$q025[i,]),
#            1:6,t(beta.trend.nat$q975[i,]), col="grey", lwd = 3)
#   
#   segments(1.1:6.1,t(beta.trend.nat.below$q025[i,]),
#            1.1:6.1,t(beta.trend.nat.below$q975[i,]), col="palegreen", lwd = 3)
#   segments(1.2:6.2,t(beta.trend.nat.below$q025[i,]),
#            1.2:6.2,t(beta.trend.nat.below$q975[i,]), col="goldenrod", lwd = 3)
#   
#   points(1:6,beta.trend.nat$mean[i,],pch= 20, col="darkgrey")
#   points(1:6,beta.trend.nat$obs[i,])
#   mtext(3, text = sub("_", " ",species[impsp[i],"tip"]), cex = 0.8)
# }
# 
# plot(1:10,1:10,type = "n", ann=F,axes=F)
# legend('left', legend = c("observed", "permute with rare occurrences", 
#                           "permute with lower abundances","permute with all plots"),
#        lty=c(NA, rep("solid",3)), lwd=c(0,5,5,5), pch=c(21, NA, NA, NA),
#        col=c("black", "grey", "palegreen", "goldenrod"), bty= "n",
#        cex = 0.8)
# 
# ## plot Standardized effect sizes
# x11()
# par(mfrow = c(4,3), mar=c(4,4,1,1))
# for (i in 1 : length(impsp)){
#   cols <- c(NA, "black") [ (beta.trend.nat$P.beta[i,] <=0.025) +1]
#   x = 1:6
#   y =(beta.trend.nat$obs[i,] - beta.trend.nat$mean[i,])/ beta.trend.nat$sd[i,]
#   y[1]=0
#   plot(x,y,ylim= c(-5,4),xlim= c(1,6.3),ann=F, las = 1, pch = 21, col ="black",bg = cols, type= "b")
#   
#   cols <- c(NA, "forestgreen") [ (beta.trend.nat.below$P.beta[i,] <=0.025) +1]
#   x = 1.1:6.1
#   y =(beta.trend.nat.below$obs[i,] - beta.trend.nat.below$mean[i,])/ beta.trend.nat.below$sd[i,]
#   par(new=T)
#   plot(x,y,ylim= c(-5,4),xlim= c(1,6.3),pch = 21, col ="forestgreen",bg = cols, type= "b",ann=F, axes=F)
#   
#   cols <- c(NA, "goldenrod") [ (beta.trend.nat.permute$P.beta[i,] <=0.025) +1]
#   x = 1.2:6.2
#   y =(beta.trend.nat.permute$obs[i,] - beta.trend.nat.permute$mean[i,])/ beta.trend.nat.permute$sd[i,]
#   par(new=T)
#   plot(x,y,ylim= c(-5,4),xlim= c(1,6.3),pch = 21, col ="goldenrod",bg = cols, type= "b",ann=F, axes = F)
#   
#   abline(h=0, lty="dotted", col="darkgrey")
#   mtext(3, text = impsp[i], cex = 0.8)
# }
# 
# plot(1:10,1:10,type = "n", ann=F,axes=F)
# legend('left', legend = c("permute with rare occurrences", 
#                           "permute with lower abundances","permute with all plots"),
#        lty=rep("solid",3),  pch=c(21, 21, 21),
#        col=c("black", "palegreen", "goldenrod"),
#        bty= "n",cex = 0.8)


spnames = rownames(glmSRnat.overall$mean)
par(mfrow = c(8,10), mar=c(0,0,2,1), oma=c(3,3,2,1))
ylim=c(-5,4)
for (i in 1 : length(spnames)){
  sp <- spnames[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ ((gamma.trend.nat$P.gamma[i,] <=0.025) | (gamma.trend.nat$P.gamma[i,] >=0.975) ) +1]
  
  # create x axis for the 5 abundance classes
  x = c(0:5)
  
  # create y axis with the gamma richness Standardized effect size
  y =(gamma.trend.nat$obs[i,] - gamma.trend.nat$mean[i,])/ gamma.trend.nat$sd[i,]
  y[1]=NA
  
  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)
  
  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
#   if (i %in% c(7:9)) text(y=-6, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)
  
#   if ( i %in% c(1,4,7)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)  
   abline(h=0,lty="dotted")
#   
  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)
  
  # Add species name
  mtext(3, text=paste(species[sp, "Genus"], "\n",species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.6, line=0, las = 1)
  
  # Y axis label
  if ( i %in% c(1,4,7)) {  
    mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
  }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=2, las = 1, outer=T)
mtext(2, text=c("Variation in native gamma richness"), adj=0.5, line=2, las = 0, outer=T)



## betas
spnames = rownames(glmSRnat.overall$mean)
par(mfrow = c(8,10), mar=c(0,0,2,1), oma=c(3,3,2,1))
ylim=c(-5,4)
for (i in 1 : length(spnames)){
  sp <- spnames[i]
  # identify significant classes in black
  cols <- c(NA, "black") [ ((geta.trend.nat$P.beta[i,] <=0.025) | (beta.trend.nat$P.beta[i,] >=0.975) ) +1]
  
  # create x axis for the 5 abundance classes
  x = c(0:5)
  
  # create y axis with the gamma richness Standardized effect size
  y =(beta.trend.nat$obs[i,] - beta.trend.nat$mean[i,])/ beta.trend.nat$sd[i,]
  y[1]=NA
  
  # plot background
  plot(x,y,ylim=ylim, xlim = c(0.5, 5.5), las = 1,type= "p", xaxt = "n", yaxt="n", ann=F,
       pch = 21, col ="black",bg = cols)
  
  axis(1,at = 1:5, label = rep("",5), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  #   if (i %in% c(7:9)) text(y=-6, x = 1:5, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)
  
  #   if ( i %in% c(1,4,7)) axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)  
  abline(h=0,lty="dotted")
  #   
  #plot points and lines
  par(new=T)
  plot(x,y, bg = cols, pch=21, type = "b",ylim=ylim, xlim = c(0.5, 5.5), xaxt = "n", yaxt="n", ann=F)
  
  # Add species name
  mtext(3, text=paste(species[sp, "Genus"], "\n",species[sp, "Species"], sep="") ,
        font = 3, outer= F,adj=0.9, cex=0.6, line=0, las = 1)
  
  # Y axis label
  if ( i %in% c(1,4,7)) {  
    mtext(2, text="SES", ,adj=0.5, cex=0.8, line=1.5, las = 0)
  }
}

mtext(1, text=c("Abundance class"), adj=0.5, line=2, las = 1, outer=T)
mtext(2, text=c("Variation in native gamma richness"), adj=0.5, line=2, las = 0, outer=T)