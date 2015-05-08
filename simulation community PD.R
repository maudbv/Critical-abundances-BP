### community simulation of faith PD bias
setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R")
library(ape)
library(picante)
load("saved Rdata/lmaphy.Rdata")

#### Example simulation

# generate 50 paired communities (25 invaded, 25 uninvaded) with different average species richness
plotdata <- data.frame( sitenames = paste("plot", 1:25, sep="_"),
                        status = c(rep ("inv", 25) , rep("non-inv", 25)),
                        SR = c( rpois(25, lambda = 9), rpois(25, lambda = 16) ))
plotdata$plotnames <- paste(plotdata$sitenames, plotdata$status, sep= "-")


### Simulated species pool of 70 species from phylogeny
# (using existing phylogeny of vascular species)
species.names <- sample(lmaphy$tip.label, 70)
phy <- drop.tip(lmaphy, tip = lmaphy$tip.label [!lmaphy$tip.label %in% species.names])



# sp1 is the invasive sp: presetn at more than 14/16 sites in inv, and absent from non inv
comm.matrix <- matrix(0, nrow= 50, ncol=70, dimnames = list( plotdata$plotnames, species.names))

# invader species is present in all invaded sites with abundance >=14
# and absent in non invaded
comm.matrix[, 1] <- c(sample(14:16, 25, replace =T), rep(0,25))

# other species present in each site have abundance from 1 to 16 :
for (l in 1:50) {
  sp <- sample(species.names[-1], plotdata$SR[l])
  for (i in sp) {
    comm.matrix[l, i] <- sample(1:16,1)
  }}

# calculate community evenness and  mean native Abundance
plotdata$GS <- diversity(comm.matrix / rowSums(comm.matrix), index = "simpson")
plotdata$mAb <- rowMeans(comm.matrix[,-1])


## Calculating phylogenetic diversity indices :
faithPD <- pd(comm.matrix, phy)

## calculate the Mean Phylogenetic distance :
modif <- comm.matrix
modif[26:50,1] <- 1   ### add the invasive species in the calculation of the MPD
MPD <- mpd(modif,cophenetic(phy) )

# Minimum distance to the invasive :
d <- cophenetic(phy) [species.names[-1],species.names[1]]
minPD <- apply(comm.matrix, 1, function(x) min(d[x[-1]>0]) )

# results:
PDresults <- data.frame(plotdata, faithPD, MPD=MPD, minPD = minPD)


### Graphs

par(mar=c(4,1,2,1), oma=c(0,0, 0, 0))
  plot(phy, cex=0.5,, show.node.label=F, type="fan")
  
  # We do have difference in SR:
  boxplot(SR ~ status, data = plotdata, ylab="Species Richness")
  
  # Faith PD is biased by SR:
  boxplot(GS ~ status , data = plotdata, ylab="Evenness") ### Non invaded plots are by definition more even in abundances!

  # PD is correlated to SR:
  plot(PD ~ SR, PDresults)
  boxplot(PD ~ status , data = PDresults, ylab="Faith's PD") #no wonder PD also decreases with invasion
  
  # MPD not correlated
  plot(MPD ~ SR, PDresults)
  boxplot(MPD ~ status , data = PDresults, ylab="MPD")
  
  # minPD not too correlated
plot(minPD ~ SR, PDresults)
quartz()  
boxplot(minPD ~ status , data = PDresults , ylab="minPD")


### Loss in phylogenetic diversity ? = PD(non-inv) - PD(inv) ?

PDresults$PDloss <-  (PDresults[1:25,"PD"] - PDresults[26:50,"PD"])/PDresults[1:25,"PD"] 
plot(PDloss ~ SR, PDresults[1:25,])




#### Simulation function
sim.random.phylo <- function(nrep = 100, paired =T,  nested = T ) {

  simul.results <- matrix(NA, nrow= nrep, ncol=7)
  colnames(simul.results ) <- c("corr.PD-SR", "corr.MPD-SR", "corr.minPD-SR",
                                "effect.PD", "effect.MPD", "effect.minPD",
                                "loss.PD-SR")

  for (k in 1:nrep) {
    
    # 25 invaded and 25 non invaded sites for 1 invasive species (ex : carpobrotus)
    plotdata <- data.frame( sitenames = paste("plot", 1:25, sep="_"),
                            status = c(rep("non-inv", 25), rep ("inv", 25) ))

    ## simulating species richness based on average SR values from previous Vila et al 2006 paper
    if(nested == F) plotdata$SR = c( rpois(25, lambda = 16), rpois(25, lambda = 9))

    ## simulating systematic nested average decrease of 25-9 = 6 species
    if(nested == T) {
      plotdata$SR = c( rpois(25, lambda = 16), rep(NA,25))
      plotdata$SR[26 : 50] = sapply(plotdata$SR[1:25], function(x){
       y <- x - rpois(1, lambda = 4)
       if (y<4) y <- min(4,x, na.rm=T)
       return(y)
    })
    }
    plotdata$plotnames <- paste(plotdata$sitenames, plotdata$status, sep= "-")


    ### Simulated species pool of 70 species from phylogeny
    # (using existing phylogeny of vascular species)

    species.names <- sample(lmaphy$tip.label, 70)
    phy <- drop.tip(lmaphy, tip = lmaphy$tip.label [!lmaphy$tip.label %in% species.names])


    # sp1 is the invasive sp: presetn at more than 14/16 sites in inv, and absent from non inv
    comm.matrix <- matrix(0, nrow= 50, ncol=70, dimnames = list( plotdata$plotnames, species.names))

    # invader species is present in all invaded sites with abundance >=14
    # and absent in non invaded
    comm.matrix[, 1] <- c(rep(0,25),sample(14:16, 25, replace =T))

    
        
    if (paired == F) {
      # species present in each site have abundance from 1 to 16 :
      for (l in 1:50) {
        sp <- sample(species.names[-1], plotdata$SR[l])
        for (i in sp) {
          comm.matrix[l, i] <- sample(1:16,1)
        }}
      
      
    }

    if (paired ==T) {
      # species present in each site have abundance from 1 to 16 :
      for (j in 1:25) {
        
        # species present in uninvaded sites
        sp <- sample(species.names[-1], plotdata$SR[j])
        for (i in sp) {
          comm.matrix[j, i] <- sample(1:16,1)
        }
        
        # species present in invaded sites
        if ( plotdata$SR[j + 25] <= plotdata$SR[j]) sp.inv <-  sample(sp, plotdata$SR[j + 25])
        if ( plotdata$SR[j + 25] > plotdata$SR[j]) {
          sp.inv <-  c(sp,sample(species.names[-1], plotdata$SR[j + 25] - plotdata$SR[j]))
        }
        
        for (ii in sp.inv) {
          comm.matrix[j + 25, ii] <- sample(1:16,1)
        }
      }
    }
    
    
    # calculate community evenness and  mean native Abundance
    plotdata$GS <- diversity(comm.matrix / rowSums(comm.matrix), index = "invsimpson")
    plotdata$mAb <- rowMeans(comm.matrix[,-1])


    ## Calculating phyloenetic diversity indices
    faithPD <- pd(comm.matrix, phy)

    modif <- comm.matrix
    modif[26:50,1] <- 1
    MPD <- mpd(modif,cophenetic(phy) )

    d <- cophenetic(phy) [species.names[-1],species.names[1]]
    minPD <- apply(comm.matrix, 1, function(x) min(d[x[-1]>0]) )


    PDresults <- data.frame(plotdata, faithPD, MPD=MPD, minPD = minPD)

    ## Loss in phylogenetic diversity ? = PD(non-inv) - PD(inv) ?
    PDresults$PDloss <-  (PDresults[1:25,"PD"] - PDresults[26:50,"PD"])/PDresults[1:25,"PD"] 
    plot(PDloss ~ SR, PDresults[1:25,])
    
    
  ## store stats
    simul.results [k, "corr.PD-SR"] <- cor.test(PDresults$PD, PDresults$SR)$est
    simul.results [k, "corr.MPD-SR"] <- cor.test(PDresults$MPD, PDresults$SR)$est
    simul.results [k, "corr.minPD-SR"] <- cor.test(PDresults$minPD, PDresults$SR)$est

    simul.results [k, "effect.PD"] <- summary(lm(PDresults$PD ~ PDresults$status))$r.squared
    simul.results [k, "effect.MPD"] <- summary(lm(PDresults$MPD ~ PDresults$status))$r.squared
    simul.results [k, "effect.minPD"] <- summary(lm(PDresults$minPD ~ PDresults$status))$r.squared

    simul.results [k, "loss.PD-SR"] <- cor.test(PDresults$PDloss[1:25], PDresults$SR[1:25])$est
     
    

    if ( k == 1) {

      plot(phy, cex=0.7)

      # We do have difference in SR :
      boxplot(SR ~ status, data = plotdata)

      # Faith PD is biased by SR :
      boxplot(GS ~ status , data = plotdata) ### Non invaded plots are by definition more even in abundances!

      plot(PD ~ SR, PDresults)
      boxplot(PD ~ status , data = PDresults) #no wonder PD also decreases with invasion

      # MPD not correlated
      plot(MPD ~ SR, PDresults)
      boxplot(MPD ~ status , data = PDresults)

      # minPD not too correlated
      boxplot(minPD ~ status , data = PDresults)
    }

  }
  return(simul.results)
}

 
 paired.sim <- sim.random.phylo(100)


Z = stack(as.data.frame(paired.sim[,1:3]))
quartz()
boxplot(values ~ ind, Z)


Z = stack(as.data.frame(paired.sim[,4:6]))
quartz()
boxplot(values ~ ind, Z)

Z = stack(as.data.frame(paired.sim[,7]))
quartz()
boxplot(values ~ ind, Z)


