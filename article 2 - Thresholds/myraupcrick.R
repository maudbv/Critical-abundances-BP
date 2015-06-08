###  Dissimilarity accounting for species richness
## null model based on lottery model similar to Chase et al. 2011 
## but coded like in Bernard-Verdier et al. 2012
##  raup crick

myraupcrick<- function(com,  nreps = 999) {
  
  bcom <- ceiling(com/10) 
bcom = bcom [rowSums(bcom >0),]

  # Community and sample sizes
  gamm<- length(colSums(bcom)>0)
  occur <- colSums(bcom)/sum(bcom) # mean relative abundance of species in the dataset
  n_sites <- nrow(bcom)
  
  # triangular matrix for formating
  tri <- matrix(FALSE, dim(bcom)[1], dim(bcom)[1])
  tri <- row(tri) > col(tri) # select the lower triangle of the matrix, excluding diagonal values

  # initiate random null model values of raup crick
  random.indices <- matrix(NA, ncol = nreps+1, nrow = sum(tri)) 
  
  # calculate nuimber of shared species for each pair of sites = Raup-crick
  random.indices[,1] <- tcrossprod(as.matrix(bcom))[tri]  
  
  # Simulate nreps random community matrices and calculate nb. of shared species 
  for (k in 1:nreps)
  {
  rcom <-  data.frame(t(apply(bcom, 1, function(x) {
      a <- sum(x)
      y <- rep(0, gamm)
      y[sample(1:gamm, a, replace = F, prob <- occur)] <- 1
      return(y)
    })))

  random.indices[,1+k] <- (tcrossprod(as.matrix(rcom)))[tri]
}


# calculate raup-crick dissimilarities as percentile of null distribution of number of shared sp.
 rc <- apply(random.indices,1, function(x) sum(x[1]<=x))
#rc <-rowSums(random.indices[,1] <= random.indices)

rm(random.indices)

# transform into distance matrix
D <- matrix(NA,  dim(bcom)[1], dim(bcom)[1])
D[tri] <- rc
# rm(rc)


Dt <- matrix(NA,  dim(bcom)[1], dim(bcom)[1])
Dt[tri] <- rc.tmp
# rm(rc)

return(rc)
}
