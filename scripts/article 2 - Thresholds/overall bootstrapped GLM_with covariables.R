### Function calculating overall bootrapped GLM tests of critical abundances

# DESCRIPTION:
# This function uses a systematic vegetation survey of species abundances (categrotical abundance classes), and quantifies the variation of a community response variable (typically, native species richness) to the increasing abundance of each focal species. 
# Focal species can be any species in the vegetation survey which matches the criteria of min.occur (i.e. only species which are widespread enough, and cover a sufficient range of abundances for the statistics to be meaningful).
# Variations of the response variable are quantified using a GLM (with poisson distribution for species richness data), including environmental covariables in addition to the focal species abundance. Significance of GLM coefficients for each abundance class are tested using a bootstrapping of the community plots performed in a previous function (informed by the input object: boot.ind)
# Coefficients of the GLM for each abundance level of the focal species are then used to detect a "Critical abundance" = the first abundance level at which coefficients are significantly negative (compared to bootstrap), and followed by only negative coefficients

# OPTIONS for this version:
# Adding flexible cofactors including BLDG distances
# Keeping all cofactors to calculate the coefficients and deviance explanations

## INPUTS:
# db :        dataframe with surveys of abundance in each plot, and columns of cofactors for each plot
# covar:      vector of column names of the environmental covariables for each plot (columns in db)
# boot.ind:   output from bootstrapping function
# variable:   column name of the response variable (in db), typically SRnat (native richness)
# min.occur:  minimum number of occurences of a focal species in a given abundance class
# min.class:  minimum number of classes a focal species needs to appear in sufficient occurrences (min.occur)
# nreps:      number of repetitions of bootstrapping

## OUTPUT is a list of the following dataframes:
# glms:         GLM statistics for each focal species 
# boot.mean:    means of response variable across bootstrapped datasets for each abundance/focal species 
# boot.sd:      standard deviations of the variable across bootstrapped datasets
# CIlow:        2.5 percentile of response variable in bootstrapped dataset
# CIhi:         97.5 percentile of response variable in bootstrapped dataset
# covar.tab:    observed coefficients and GLM statistics for each covariable in each GLM (one per focal species)
# crit.vals:    Critical abundance values detected for each focal species
# spearman:     statistics of spearman's sign rank test for each focal species 
# n.obs:        number of observed occurrences for each abundance class/focal species
# mean.values:  means of the response variable across observed plots for each abundance class/focal species
# dif = dif:    difference between mean response variable between consecutive abundance classes for each species
# est= est:     observed coefficients of GLMs for each abundance class/focal species.
# z=z:          standardized effect size for each difference of means
# P= P:         P value of the GLM coefficients for each abundance class/focal species.


# This function was written by Maud Bernard-Verdier - updated March 2018
# ______________________________________________________________________________________________________


glm.overallboot<- function(db=databp[databp$PlotName %in% unimprovedgrasslands,], covar = c( 'DEM_10', "SLOPE",'Northern', 'SRali'),
                           boot.ind  = boot.indices , 
                           variable='SRnat', min.occur =5,  min.class = 2, nreps = 9) {
  
  # Parallel computation when possible:
  library(doParallel)
  cl<-makeCluster(6)
  registerDoParallel(cl)

  # check that we have covariables:
  stopifnot(length(covar) >= 1,  is.character(covar) )
  
  ### Check that boot.output is the right size:
  stopifnot(min.occur == boot.ind$min.occur,
            min.class == boot.ind$min.class,
            dim(db)[1] == boot.ind$db.size)
  
  
  ## Selecting focal species verifying conditions :
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=min.class)
                    &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
  db.modif <- db[which(db$SpeciesCode %in% a),]
  
  # list of focal species:
  sp.names <- unique(db.modif$SpeciesCode)
  
  
  #### Calculating GLM for observed and bootstrap datasets for each species
  
  # initiate result dataframes :
  crit.vals <- crit.vals.P  <- crit.vals.CI <-  data.frame(matrix(NA,
                                                                  nrow=length(sp.names),
                                                                  ncol= nreps +1,
                                                                  dimnames=list(sp.names,1:(nreps+1))))
  
  glms <- data.frame(
    matrix(NA, nrow=length(sp.names), ncol= 7 + (length(covar)*2), 
           dimnames=list(sp.names,
                         c( "df", "resid.dev","null.dev", "dev.ratio","aic.null",
                            paste("aic", c(covar, "abun"), sep = "."),
                            paste("dev", c(covar, "abun"), sep = ".")
                            )
                         )
           )
    )
  
  
  init <-  data.frame(matrix(NA,
                             nrow=length(sp.names),
                             ncol= 5,
                             dimnames=list(sp.names,c("c2", "c3", "c4", "c5", "c6"))))
  
  dif <-   est<- P <- z <- init
  
  n.obs <- mean.values <- data.frame(matrix(NA,
                                            nrow=length(sp.names),
                                            ncol= 6,
                                            dimnames=list(sp.names,c("C1","c2", "c3", "c4", "c5", "c6"))))
  
  spear <-data.frame(matrix(NA,
                            nrow=length(sp.names),
                            ncol= 2,
                            dimnames=list(sp.names, c("rho" ,"p.val"))))
  
  impact.size <-  data.frame(matrix(NA,
                                    nrow=length(sp.names),
                                    ncol= 9,
                                    dimnames=list(sp.names,
                                                  c("th","th.CI","prevalence", "nb.plot.impact","prop.plot.impact","mean.dif","wtd.mean.dif", "th.dif", "max.dif"))))
  

  bootCI.low <-bootCI.hi <- boot.mean <- boot.sd <-  init
  
  
  covar.tab <- list()
  for (k in 1:length(covar)) covar.tab[[k]] <- data.frame(matrix(NA, nrow=length(sp.names), ncol= 6, dimnames=list(sp.names,c("coef.glm", "P.coef", "mean", "sd", "lowCI", "hiCI"))))
  names(covar.tab) = covar
  
  
  ## Calculations looping on all focal species:

  for (i in 1:length(sp.names) ) {
    
    sp <- sp.names[i]  # select species name
    print(paste(i, ":", sp, "(",Sys.time(),")"))
    
    ### FIRST : calculate GLMs and other stats for observed dataset :
    sp.dat <- db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",variable,'PlotName', covar)] # select occurrences of the species
    names(sp.dat) <- c('abun','var','PlotName',covar)[1:dim(sp.dat)[2]]
    
    abun <- sort(as.numeric(as.character(na.omit(unique(sp.dat$abun))))) ## list of abundance classes for species i
    n.obs[i,] <- sapply(1:6, FUN=function(j) length(sp.dat[which(sp.dat$abun==j),"var"]))
    
    #calculate differences in response variables between consecutive abundance classes:
    dif[i,] <-  sapply(2:6, FUN=function(j) mean(sp.dat[which(sp.dat$abun==j),"var"], na.rm=T) - mean (sp.dat[which(sp.dat$abun==min(abun)),"var"], na.rm=T))
    mean.values[i,] <- sapply(1:6, FUN=function(j) mean(sp.dat[which(sp.dat$abun==j),"var"], na.rm=T))
    
    # spearmnan test across all classes
    s <- cor.test(sp.dat$var,sp.dat$abun, method="spearman")
    spear[i,] <-  c(s$estimate, s$p.value)
    
    # GLM test with AIC 
    fits <- as.list(1:(length(covar) + 2))
    names(fits) <- c("f0", paste("f",1:length(covar), sep = ""), "ff")
    
    # empty GLM (no explanatory variable)
    fits$f0 <-  glm(as.formula(paste("var ~ 1")), data = sp.dat,  family=poisson(log))
    
    # Successive GLMs with increasing number of covariables:
    for (k in 1:length(covar)) {
      fits[[1+k]] <-  glm(as.formula(paste("var ~", paste(covar[1:k], collapse = " + "))),
                          data = sp.dat,  family=poisson(log))
    }  
    
    # Final GLM with all covariables and the focal species abundance (aka: what we want to test)
    fits$ff <-  glm(as.formula(paste("var ~", paste(covar[1:length(covar)], collapse = " + "),"+ as.factor(abun)")), data = sp.dat,  family=poisson(log))
      
    # Quantifying co-variable contributions to AIC reduction
    aics <- sapply(1:length(names(fits)), function(x) AIC(fits[[x]]))
    names(aics) <- names(fits)

    # store full model:
    f <- fits$ff
    
    # Storing GLM results:
    glms[i,] <-  c(df= f$df.resid,
                   resid.dev= f$dev, 
                   null.dev =f$null.deviance, 
                   dev.ratio= (f$null.deviance - f$dev)/f$null.deviance, 
                   aics,
                   anova(f)$Deviance[-1])
    
    #extract the observed coefs for species i:
    obs.coef <- rep(NA,5)
    for (j  in abun[-1])  {
      obs.coef[j-1]  <- summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]
      z[i,j-1]  <- summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )),3]
      P[i,j-1]  <- summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )),4]
    }
    names(obs.coef) <- c("c2", "c3", "c4", "c5", "c6")
    est[i,] <- obs.coef
    
    # observed GLM coefficients and stats for all the covariables
      for ( k in c(covar)) {
        covar.tab[[k]]$coef.glm [i] <- summary(f)$coef [grep(k, rownames(summary(f)$coef )), 1]
        covar.tab[[k]]$P.coef [i] <- summary(f)$coef [grep(k, rownames(summary(f)$coef )), 4]
      }
    
  
    ## SECOND : recalculate coefs for each of the *nreps* bootstrapped datasets:
  
    # bootstrapping of coefficients nreps times
    boot.coef <-  foreach(r=1:nreps, .combine='rbind') %dopar%  {
      
      ## extract new datasets from original dataset using the bootstrapped line numbers:
      boot.db <- db.modif [boot.ind$index[boot.ind$index$nrep == r, "ind" ], c("SpeciesCode", "abun",variable,'PlotName', covar)]
      
      # print(paste(r, ":","(",Sys.time(),")"))
      dat <- boot.db[as.character( boot.db$SpeciesCode)==sp, ] # select occurrences of the species
      names(dat) <- c("sp",'abun','var','PlotName', covar)[1:dim(dat)[2]]
      rm(boot.db)
      abun <- sort(as.numeric(as.character(na.omit(unique(dat$abun))))) ## list of abundance classes for species i
      
      covar.coefs <- rep(NA, length(covar))
      names(covar.coefs) = covar
      
      ### Calculate GLM model if possible (sufficient abundance classes)
      if ("1" %in% abun & length(abun) >1 )
      {
          f <-  glm(as.formula(paste("var ~",paste(covar, collapse = " + ")," + as.factor(abun)")), data = dat,  family=poisson(log))
          ### store coefficients for each bootstrapped sample k
          for (k in covar) {
            covar.coefs[k] <- summary(f)$coef [grep(k, rownames(summary(f)$coef )), 1]
          }
        
        coef.string <- rep(NA,5)
        for (j  in abun[-1]) {
          if (length(grep(paste(j), rownames(summary(f)$coef )))>0){
            coef.string[j-1] <- summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]
          }
        }
        
        # if (all(is.na(coef.string))) coef.string <- rep(NA,5)
        
      }
      
      return(c(covar.coefs, coef.string))
    }
    
    # assemble observed and bootstrapped coefs
    coefs <- rbind(obs = c(sapply(covar,FUN=function(k) covar.tab[[k]]$coef.glm [i]),obs.coef),
                   boot.coef)
    
    
    
    
    #### THIRD: statistics for species i

      # Statistics for covariables:
      for ( k in covar) {
        covar.tab[[k]]$mean[i] <- mean(coefs[,k], na.rm = T)
        covar.tab[[k]]$sd[i] <- sd(coefs[,k], na.rm = T)
        covar.tab[[k]]$lowCI[i] <-quantile(coefs[,k] , probs = 0.025,na.rm =T)
        covar.tab[[k]]$hiCI[i] <-quantile(coefs[,k] , probs = 0.975,na.rm =T)
      }
    
    # detect minimal critical value (simple start of negative trend)
    coefs <- coefs[,-which(colnames(coefs) %in% covar)]
    
    crit.vals[i,] <- sapply(1:(nreps+1), function(x, n.obs) {
      crit <- NA
      co <-coefs [x,]
      
      if (!all (is.na(co))) {
        ab = which(!is.na(co))+1  
        neg <-   which(co<0 ) +1
        
        candidates <- which(co<0 & (n.obs[2:6] >= min.occur) ) +1   ## good candidates: negative coefs, sufficient number of obs
        
        if (length(candidates)>=1) {
          y <- candidates [sapply(candidates, FUN= function(l) {
            c1 <- ( if ( l+1 <= max(abun)) all(((l+1):max(ab)) %in% neg) # all higher classes have negative diferences
                    else c1 =F)
            return(c1)
          })]
          if (length(y)>=1) crit <- min( y, na.rm=T)
        }
      }
      return(crit)
    }, n.obs  = n.obs[i,])
    
    
    # calculate CI for bootstrapped coefs
    boot.mean [i,] <-apply(coefs, 2, mean,na.rm =T)
    boot.sd [i,] <-apply(coefs, 2, sd,na.rm =T)
    
    bootCI.low [i,] <-apply(coefs, 2, quantile, probs = 0.025,na.rm =T)
    bootCI.hi [i,] <-apply(coefs, 2, quantile, probs = 0.975,na.rm =T)
    
  }
  
   stopCluster(cl) # Only if parallel computing
  
  return(list(glms=glms,boot.mean = boot.mean, boot.sd=boot.sd, CIlow = bootCI.low, CIhi = bootCI.hi,
              covar.tab = covar.tab,
              crit.vals = crit.vals, spearman=spear, n.obs = n.obs,
              mean.values = mean.values,
              dif = dif,est= est, z=z,P= P))
  
}