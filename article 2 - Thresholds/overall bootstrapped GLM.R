### overall bootrapped GLM test : Elena's advice


###-----------Creating bootstrap samples of the dataset :
bootstrap.dataset<- function(db=databp[databp$PlotName %in% realgrasslands,], 
                             min.occur =5,  min.class = 2, nreps = 999) {
  
  ######### selecting species verifying conditions : 
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=min.class) 
                    &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
  db.modif <- db[which(db$SpeciesCode %in% a),] 
  
  # list of species =to be targeted in analysis :
  sp.names <- unique(db.modif$SpeciesCode)
  
  
  ######### Bootstrapping the dataset : resampling plots with replacement
  
  # Unique set of plots to be resampled 
  d = as.character(unique( db.modif$PlotName))
  
  # resampling the list of plot names with replacement *nreps* times
  boots <- sapply(1:nreps,FUN = function (k) sample(d, replace=TRUE))
  
  return(list(boots =boots,db.size=dim(db)[1],min.occur =min.occur,  min.class = min.class, nreps = nreps))
}

extract.indices <- function(boot.output, db = db) {
  
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=boot.output$min.occur)>=boot.output$min.class) 
                    &  table(db$SpeciesCode, db$abun)[,1]>=boot.output$min.occur))
  db.modif <- db[which(db$SpeciesCode %in% a),] 
  
  indices.table=NULL
  
  for (r in 1:boot.output$nreps){
    ## extract new datasets from original dataset using the bootstrapped line numbers:
    ind <-unlist(lapply(1:length(boot.output$boots[,r]), function (k)  which(db.modif$PlotName == boot.output$boots[k,r])))
    nrep = rep(r, length(ind))
    indices.table =rbind(indices.table, cbind(nrep, ind))
    print(paste(r, ":","(",Sys.time(),")"))
  }
  
  return(list(index = as.data.frame(indices.table), nreps = boot.output$nreps, min.class = boot.output$min.class, min.occur = boot.output$min.occur, db.size = boot.output$db.size))
}


###-----------Calculating GLMs, critical values and impact size for bootstrapped samples

glm.overallboot<- function(db=databp[databp$PlotName %in% realgrasslands,], 
                           boot.indices  = boot.indices , sp.target = NA,
                           variable='SR', min.occur =5,  min.class = 1, nreps = 999) {
  
  library(doParallel)
  cl<-makeCluster(6)
  registerDoParallel(cl)
  
  ### Check that boot.output is the right one
  stopifnot(min.occur == boot.indices$min.occur,
            min.class == boot.indices$min.class,
            dim(db)[1] == boot.indices$db.size)
  
  
  ######### selecting species verifying conditions : 
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=min.class) 
                    &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
  db.modif <- db[which(db$SpeciesCode %in% a),] 
  
  # list of species =to be targeted in analysis :
  sp.names <- unique(db.modif$SpeciesCode)
  
  
  #### Calculating GLM for observed and bootstrap datasets for each species
  
  # initiate result dataframes :
  crit.vals <- crit.vals.P  <- crit.vals.CI <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= nreps +1,
                                                                  dimnames=list(sp.names,1:(nreps+1))))
  # coefs.C2 <- coefs.C3 <- coefs.C4 <- coefs.C5 <- coefs.C6 <- crit.vals
  # boot.coefs <- list(coefs.C2,coefs.C3,coefs.C4,coefs.C5,coefs.C6)
  
  
  
  init <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 5,
                             dimnames=list(sp.names,c("c2", "c3", "c4", "c5", "c6"))))
  dif <-   est<- P <- z <- init
  
  n.obs <- mean.values <- data.frame(matrix(NA, nrow=length(sp.names), ncol= 6,
                                            dimnames=list(sp.names,c("C1","c2", "c3", "c4", "c5", "c6"))))
  
  glms <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 3,
                           dimnames=list(sp.names, c( "df", "resid.dev","dev.ratio"))))
  
  spear <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 2,
                            dimnames=list(sp.names, c("rho" ,"p.val"))))
  
  impact.size <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 9,
                                    dimnames=list(sp.names, c("th","th.CI","prevalence", "nb.plot.impact","prop.plot.impact",
                                                              "mean.dif","wtd.mean.dif", "th.dif", "max.dif"))))
  
  impact.spread <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 5,
                                      dimnames=list(sp.names, c("th","th.CI","prevalence", "nb.plot.impact","prop.plot.impact"))))
  
  
  bootCI.low <-bootCI.hi <- boot.mean <- boot.sd <-  init 
  
  
  
  # looping on species
  for (i in 1:length(sp.names) ) {
    
    sp <- sp.names[i]  # select species name
    print(paste(i, ":", sp, "(",Sys.time(),")"))
    
    ### FIRST : calculate GLM for observed dataset :
    sp.dat <- db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",variable,'PlotName')] # select occurrences of the species
    names(sp.dat) <- c('abun','var','PlotName')[1:dim(sp.dat)[2]]
    
    abun <- sort(as.numeric(as.character(na.omit(unique(sp.dat$abun))))) ## list of abundance classes for species i
    n.obs[i,] <- sapply(1:6, FUN=function(j) length(sp.dat[which(sp.dat$abun==j),"var"]))
    
    #calculate diference in class mean SR
    dif[i,] <-  sapply(2:6, FUN=function(j) mean(sp.dat[which(sp.dat$abun==j),"var"], na.rm=T) - mean (sp.dat[which(sp.dat$abun==min(abun)),"var"], na.rm=T))
    mean.values[i,] <- sapply(1:6, FUN=function(j) mean(sp.dat[which(sp.dat$abun==j),"var"], na.rm=T))
    
    # spearmnan test across all classes
    s <- cor.test(sp.dat$var,sp.dat$abun, method="spearman")  
    spear[i,] <-  c(s$estimate, s$p.value)
    
    # GLM test
    f <-  glm(sp.dat$var ~ as.factor(sp.dat$abun), family=poisson(log))
    glms[i,] <-  c(df= f$df.resid, resid.dev= f$dev,dev.ratio= (f$null.deviance -f$dev)/f$null.deviance )
    n <-  1:(length(abun)-1)
    est[i,n] <- summary(f)$coef [-1, 1]
    
    z[i,n] <- summary(f)$coef [-1, 3]
    P[i,n] <- summary(f)$coef [-1, 4]
    
    
    ### if species has a negative coefficient :
    # sp.sig <-  which(!apply( !(est<0 & n.obs[,2:6] >=5), 1, all, na.rm=T))
    ### or is in a list of targeted species :sp.target
#     if (    !all(!( est[i,]<0 & n.obs[i,2:6] >=5), na.rm=T) | sp %in% sp.target  ) {
      
#       coefs <- data.frame(matrix(NA, nrow=nreps+1, ncol= 5,
#                                  dimnames=list(1:(nreps+1) ,c("c2", "c3", "c4", "c5", "c6"))))
#       pvals <- data.frame(matrix(NA, nrow=nreps+1, ncol= 5,
#                                  dimnames=list(1:(nreps+1) ,c("c2", "c3", "c4", "c5", "c6"))))
#       
#       for(j  in abun[-1]) coefs[1, j-1] = summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]
#       for(j  in abun[-1]) pvals[1, j-1] = summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 4]
      
      #store coefs
      #   for(j  in abun[-1]) boot.coefs[[j-1]][1, k] = summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]
      #   

obs.coef <- rep(NA,5)
for(j  in abun[-1]) coefs[1, j-1] <- summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]
names(obs.coef) <- c("c2", "c3", "c4", "c5", "c6")

## SECOND : recalculate for each of the *nreps* bootstrapped datasets:
boot.coef <-  foreach(i=1:nreps, .combine='rbind') %dopar%  {
        
  ## extract new datasets from original dataset using the bootstrapped line numbers:
              boot.db <- db.modif [boot.indices$index[boot.indices$index$nrep == r, "ind" ], c("SpeciesCode", "abun",variable,'PlotName')]
        
        # print(paste(r, ":","(",Sys.time(),")"))
        dat <- boot.db[as.character( boot.db$SpeciesCode)==sp, ] # select occurrences of the species
        names(dat) <- c("sp",'abun','var','PlotName')[1:dim(dat)[2]]
        rm(boot.db)
        abun <- sort(as.numeric(as.character(na.omit(unique(dat$abun))))) ## list of abundance classes for species i
        
        ### Calculate GLM model if possible (sufficient abundance classes)
        if ("1" %in% abun & length(abun) >1 ) 
        {
          f <-  glm(dat$var ~ as.factor(dat$abun), family=poisson(log))
          
          # ## store coefficients for each bootstrapped sample k
          coef.string <- rep(NA,5)
          for(j  in abun[-1]) coefs[j-1] <- summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]
          print(paste(i, ":",r, ":", "(",Sys.time(),")"))
      }
      

      #### THIRD: statistics for species i 
      
      # detect minimal critical value (simple start of negative trend)
      
      crit.vals[i,] <- sapply(1:(nreps+1), function(x, n.obs) {
        crit <- NA
        co <-coefs [x,]
        ab = which(!is.na(co))+1
        neg <-   which(co<0 ) +1
        
        candidates <- which(co<0 & n.obs[2:6] >= min.occur ) +1
      
        if (length(candidates)>=1) {
          y <- candidates [sapply(candidates, FUN= function(l) {
            c1 <- ( if ( l+1 <= max(abun)) all(((l+1):max(ab)) %in% neg) # all higher classes have negative diferences
                    else c1 =F)
            return(c1)
          })]
          if (length(y)>=1) crit <- min( y, na.rm=T)
        }
        return(crit)
      }, n.obs  = n.obs[i,])
      
      
      # detect minimal critical value using GLM Pvalues
      # crit.vals.P[i,] <- sapply(1:(nreps+1), function(x, n.obs) {
      #   crit <- NA
      #   co <-coefs [x,]
      #   ps <- pvals [x,]
      #   ab = which(!is.na(co))+1
      #   neg <-   which(co<0) +1
      #   sig <-  which(ps<=0.05 & n.obs[2:6] >= min.occur) +1
      #   sig= sig [sig %in% neg]
      # 
      #   if (length(sig)>=1) {
      #     y <- sig [sapply(sig, FUN= function(l) {
      #       c1 <- ( if ( l+1 <= max(abun)) all(((l+1):max(ab)) %in% neg) # all higher classes have negative diferences
      #               else c1 =F)
      #       return(c1)
      #     })]
      #     if (length(y)>=1) crit <- min( y, na.rm=T)
      #   }
      #   return(crit)
      # }, n.obs  = n.obs[i,])
      
      # calculate CI for bootstrapped coefs
      boot.mean [i,] <-apply(coefs, 2, mean,na.rm =T) 
      boot.sd [i,] <-apply(coefs, 2, sd,na.rm =T) 
      
      bootCI.low [i,] <-apply(coefs, 2, quantile, probs = 0.025,na.rm =T) 
      bootCI.hi [i,] <-apply(coefs, 2, quantile, probs = 0.975,na.rm =T) 
      
      testCI = sign(as.numeric(bootCI.low [i,])) * sign(as.numeric(bootCI.hi [i,]))
      
      # detect minimal critical value using bootstrap CI
      crit.vals.CI[i,] <- sapply(1:(nreps+1), function(x, n.obs) {
        crit <- NA
        co <-coefs [x,]
        ps <- pvals [x,]
        ab = which(!is.na(co))+1
        neg <-   which(co<0) +1
        
        # select coefficients whose CI are above or below zero :
        sig <-  which(testCI ==T & n.obs[,2:6] >= min.occur) + 1
        
        sig= sig [sig %in% neg]
        
        if (length(sig)>=1) {
          y <- sig [sapply(sig, FUN= function(l) {
            test <- F
            if (l == max(abun, na.rm=T)) test <- T
            if (l < max(abun, na.rm=T)) test<- all(((l+1):max(ab)) %in% neg) # all higher classes have negative diferences
            return(test)
          })]
          if (length(y)>=1) crit <- min( y, na.rm=T)
        }
        return(crit)
      }, n.obs  = n.obs[i,])
      
      
      
    }
    
#   }
  
stopCluster(cl)

  return(list(glms=glms, impact.spread = impact.spread,
              boot.mean = boot.mean, boot.sd=boot.sd, CIlow = bootCI.low, CIhi = bootCI.hi,
              crit.vals = crit.vals,crit.vals.P = crit.vals.P, spearman=spear, n.obs = n.obs,
              mean.values = mean.values,
              dif = dif,est= est, z=z,P= P))
  }
}