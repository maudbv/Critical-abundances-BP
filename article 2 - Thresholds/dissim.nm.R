###  Dissimilarity accounting for species richness
## null model based on lottery model similar to Chase et al. 2011 
## but coded like in Bernard-Verdier et al. 2012
##  raup crick

dissim.nm<- function(com, nullmodel = T, nreps = 999, index = c("rc", "sm", "euc")) {
  
  
  source(file ="script/functions/SMsim.R")
  
  #transform to presence/absence
  bcom <- ceiling(com>0) 
# bcom = bcom [rowSums(bcom >0),]

  # Community and sample sizes
  gamm<- length(colSums(bcom)>0)
  occur <- colSums(bcom)/sum(bcom) # mean relative abundance of species in the dataset
  n_sites <- nrow(bcom)
  
  # triangular matrix for formating
  tri <- matrix(FALSE, dim(bcom)[1], dim(bcom)[1])
  tri <- row(tri) > col(tri) # select the lower triangle of the matrix, excluding diagonal values

  # initiate random null model values of raup crick
  random.indices <- matrix(NA, ncol = nreps+1, nrow = sum(tri)) 
  
if (index == "rc") {# calculate nuimber of shared species for each pair of sites = Raup-crick
  random.indices[,1]<- tcrossprod(as.matrix(bcom))[tri]  
}

if (index == "euc") {  #  euclidean on presence /absence
  random.indices[,1] <- as.matrix(vegdist(as.matrix(bcom), method = "euclid"))[tri]  
}

if (index == "sm") {  #  Species matching similarity  = (a + d) / (a + d + b + c)
  random.indices[,1] <-  as.matrix(1-SMsim(bcom))[tri]  
}


  # Simulate nreps random community matrices and calculate nb. of shared species 
  for (k in 1:nreps)
  {

# trying will's method : doesnt seem to work well : accounts for SR in plots, but not sp occurrence prob
#     probs <- matrix(occur, n_sites, gamm, byrow = T)   
#     sitep <- matrix(sr/gamm, n_sites, gamm, byrow = F)
#     
#     ps <- sitep * probs
#     rcom2 <- t(apply(ps,1, function(x) rbinom(gamm,1,x)))
          
  rcom <-  data.frame(t(apply(bcom, 1, function(x) {
      a <- sum(x)
      y <- rep(0, gamm)
      y[sample(1:gamm, a, replace = F, prob <- occur)] <- 1
      return(y)
    })))
  
  if (index == "rc") {# calculate nuimber of shared species for each pair of sites = Raup-crick
    random.indices[,1+k]<- tcrossprod(as.matrix(rcom))[tri]  
  }
  
  if (index == "euc") {  #  euclidean on presence /absence
    random.indices[,1+k] <- as.matrix(vegdist(as.matrix(rcom), method = "euclid"))[tri]  
  }
  
  if (index == "sm") {  #  Species matching similarity  = (a + d) / (a + d + b + c)
    random.indices[,1+k] <-  as.matrix(1-SMsim(rcom))[tri]  
  }

}


# calculate index as percentile of null distribution 
p.dist <- apply(random.indices,1, function(x) {
 if (!all(unique(x) == 0)) ind <- sum(x[1]>=x)
if (all(unique(x) == 0)) ind<- 0
return(ind)
})/(nreps+1)

# calculate effect size : makes negative values
#index <- apply(random.indices,1, function(x) (x[1] - mean(x))/sd(x))


#rc <-rowSums(random.indices[,1] <= random.indices)

rm(random.indices)

# transform into distance matrix
D <- matrix(NA,  dim(bcom)[1], dim(bcom)[1])
D[tri] <- p.dist
rm(p.dist)

return(as.dist(D))
}
