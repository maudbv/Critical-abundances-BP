    ## homemade bootstrapping function using document by John Fox 2002 Bootstrapping Regression Models
    ## http://cran.r-project.org/doc/contrib/Fox-Companion/appendix-bootstrapping.pdf
    myboot <- function (d =sp.dat, fn = f.glm, R = R, keep.size = TRUE, inextenso = F) {
     
     if (inextenso == F) {
      boot.coeff <- as.data.frame(sapply(1:R,FUN = function (k, keep.size = keep.size) {
         if (keep.size) xb <- unlist(tapply(1:length(d$abun), INDEX = d$abun, FUN = sample, replace=TRUE )) ## keeps sample size in each classs of abundance 
         if (!keep.size) xb <- sample(1:length(d$abun), replace=TRUE ) ## changes the sample size in each class of abundance
          fb <- f.glm(d = sp.dat, w = xb)
          return(fb)
        }, keep.size=keep.size))
       
        boot.coeff$t0 <- f.glm(sp.dat, rownames(sp.dat))
        boot.CI <- cbind(t0 =    boot.coeff$t0,
                         as.data.frame(t(apply(boot.coeff, 1,quantile, probs=c(0.025, 0.975)))))
        boot.CI$Pnegative <- apply(boot.coeff, 1, FUN= function(k) sum(k > 0)/(R+1))
        bias <-  rowMeans(boot.coeff) - boot.coeff$t0 
        boot.se <- sqrt(rowSums(apply(boot.coeff,2, FUN= function(k) (k -  rowMeans(boot.coeff) )^2)) / (R) )
        boot.out <- data.frame(boot.CI, bias = bias, se = boot.se ) 
        return(boot.out)
     }
     
     if (inextenso == T) {
     indices <- as.data.frame(sapply(1:R,FUN = function (k, keep.size = keep.size) {
         if (keep.size) xb <- unlist(tapply(1:length(d$abun), INDEX = d$abun, FUN = sample, replace=TRUE )) ## keeps sample size in each classs of abundance 
         if (!keep.size) xb <- sample(1:length(d$abun), replace=TRUE ) ## changes the sample size in each class of abundance
         return(xb)
       }, keep.size=keep.size))
      return(indices)
      }
  }
    
    
    boot.results <- myboot(sp.dat, fn = f.glm, R = R, keep.size = F, inextenso = T)
    
    
    