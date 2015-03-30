    ## homemade bootstrapping function using document by John Fox 2002 Bootstrapping Regression Models
    ## http://cran.r-project.org/doc/contrib/Fox-Companion/appendix-bootstrapping.pdf

   overall<- function (d = as.character(unique(db$PlotName)), R = 99) {
     boots <- sapply(1:R,FUN = function (k) sample(d, replace=TRUE)) ## resample with replacement
     indices <- cbind(original = d, boots)
      }
  
    
    
    boot.results <- myboot(sp.dat, fn = f.glm, R = R, keep.size = F, inextenso = T)
    
    
#     ### include glm calculations on each boot sample
#     
#     if (justtheindices == F) {
#       boot.coeff <- sapply(1:R,FUN = function (k, keep.size = keep.size) {
#         if (keep.size) xb <- unlist(tapply(1:length(d$abun), INDEX = d$abun, FUN = sample, replace=TRUE )) ## keeps sample size in each classs of abundance 
#         if (!keep.size) xb <- sample(1:length(d$abun), replace=TRUE ) ## changes the sample size in each class of abundance
#         fb <- f.glm(d = sp.dat, w = xb)
#         names(fb) <- sapply( names(fb), FUN= function(nc) substr(nc, nchar(nc),nchar(nc)) )
#         names(fb)[1] <- "intercept"
#         return(fb)
#       },
#       keep.size=keep.size)
#       
#       boot.coeff$t0 <- f.glm(sp.dat, rownames(sp.dat))
#       names( boot.coeff$t0) <- sapply( names( boot.coeff$t0), FUN= function(nc) substr(nc, nchar(nc),nchar(nc)) )
#       names( boot.coeff$t0)[1] <- "intercept"
#       
#       boot.coeff <-  sapply(names(boot.coeff$t0), FUN = function(nc) {
#         print(nc)
#         val <- as.numeric(lapply(boot.coeff, FUN = function(m) m[grep(nc , names(m))]))
#         return(val)
#       })
#       
#       
#       boot.CI <- t(rbind(t0 =    boot.coeff[1,],
#                          apply(boot.coeff,2, quantile, probs=c(0.025, 0.975), na.rm =T)))
#       boot.CI$Pnegative <- apply(boot.coeff,2, FUN= function(k) sum(k > 0, na.rm=T)/sum(!is.na(k)))
#       bias <-  colMeans(boot.coeff, na.rm=T) - boot.coeff[1,]
#       boot.se <- sqrt(colSums(apply(boot.coeff,2, FUN= function(k) (k -  colMeans(boot.coeff, na.rm=T)^2) / sum(!is.na(k)) )
#                               boot.out <- data.frame(boot.CI, bias = bias, se = boot.se ) 
#                               return(boot.out)
#     }
#     
#     