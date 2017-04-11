SMsim <- function(com){
  com=ceiling(com>0)
  a <- tcrossprod(as.matrix(com))  # joint presences
  d <- tcrossprod(ceiling(as.matrix(com==0))) # joint absences
  
  sr <- rowSums(com)
  sri = matrix(sr, nrow = length(sr), ncol = length(sr), byrow=F)  
  srj = matrix(sr, nrow = length(sr), ncol = length(sr), byrow=T)  
  b <- as.dist(srj - a)   # species only in site i
  c <- as.dist(sri - a)   # species only in site j
  a<- as.dist(a)
  d <- as.dist(d)
  Sim <-  (a+d) / (a +b+c +d)
  
  return(Sim)
}


### shows how SM index is influenced by the difference in sp richness => need to standardize by diff in species richness
# 
# sr <- rowSums(bcom)
# sri = matrix(sr, nrow = length(sr), ncol = length(sr), byrow=F)  
# srj = matrix(sr, nrow = length(sr), ncol = length(sr), byrow=T) 
# 
# plot(as.dist(2*abs(sri - srj)/sri + srj) , SMsim(bcom))   
# plot(as.dist(2*abs(sri - srj)/sri + srj) , vegdist(bcom, method = "euclid"))  
# plot(as.dist(2*abs(sri - srj)/sri + srj) , vegdist(bcom, method = "bray"))
# plot(as.dist(2*abs(sri - srj)/sri + srj) , vegdist(bcom, method = "jaccard"))
