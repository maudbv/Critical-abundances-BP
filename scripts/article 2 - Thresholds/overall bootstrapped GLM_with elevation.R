### overall bootrapped GLM test : Elena's advice


###-----------Calculating GLMs, critical values and impact size for bootstrapped samples + elevation as co-variable

glm.overallboot<- function(db=databp[databp$PlotName %in% unimprovedgrasslands,], covar = 'DEM_10',
                           boot.ind  = boot.indices , sp.target = NA,
                           variable='SR', min.occur =5,  min.class = 2, nreps = 9) {

  library(doParallel)
  cl<-makeCluster(6)
  registerDoParallel(cl)

  ### Check that boot.output is the right one
  stopifnot(min.occur == boot.ind$min.occur,
            min.class == boot.ind$min.class,
            dim(db)[1] == boot.ind$db.size)


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

  covar.tab <- data.frame(matrix(NA, nrow=length(sp.names), ncol= 6,
                         dimnames=list(sp.names,c("coef.glm", "P.coef", "mean", "sd", "lowCI", "hiCI"))))
  

  # looping on species
  for (i in 1:length(sp.names) ) {

    sp <- sp.names[i]  # select species name
    print(paste(i, ":", sp, "(",Sys.time(),")"))

    ### FIRST : calculate GLM for observed dataset :
    sp.dat <- db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",variable,'PlotName', covar)] # select occurrences of the species
    names(sp.dat) <- c('abun','var','PlotName','covar')[1:dim(sp.dat)[2]]
      
    abun <- sort(as.numeric(as.character(na.omit(unique(sp.dat$abun))))) ## list of abundance classes for species i
    n.obs[i,] <- sapply(1:6, FUN=function(j) length(sp.dat[which(sp.dat$abun==j),"var"]))

    #calculate diference in class mean SR
    dif[i,] <-  sapply(2:6, FUN=function(j) mean(sp.dat[which(sp.dat$abun==j),"var"], na.rm=T) - mean (sp.dat[which(sp.dat$abun==min(abun)),"var"], na.rm=T))
    mean.values[i,] <- sapply(1:6, FUN=function(j) mean(sp.dat[which(sp.dat$abun==j),"var"], na.rm=T))

    # spearmnan test across all classes
    s <- cor.test(sp.dat$var,sp.dat$abun, method="spearman")
    spear[i,] <-  c(s$estimate, s$p.value)

    # GLM test
    f <-  glm(sp.dat$var ~ as.numeric(sp.dat$covar) + as.factor(sp.dat$abun), family=poisson(log))
    glms[i,] <-  c(df= f$df.resid, resid.dev= f$dev, dev.ratio= (f$null.deviance -f$dev)/f$null.deviance )
    n <- (abun-1)[-1]
    est[i,n] <- summary(f)$coef [-(1:2), 1]
    z[i,n] <- summary(f)$coef [-(1:2), 3]
    P[i,n] <- summary(f)$coef [-(1:2), 4]


    ## SECOND : recalculate coefs for each of the *nreps* bootstrapped datasets:

    #extract the observed coefs for species i
    obs.coef <- rep(NA,5)
    for (j  in abun[-1])  obs.coef[j-1]  <- summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]
    names(obs.coef) <- c("c2", "c3", "c4", "c5", "c6")

    # observed GLM values for the covariable
    covar.tab$coef.glm [i] <- summary(f)$coef [2, 1]
    covar.tab$P.coef [i] <- summary(f)$coef [2,4]
    
    # bootstrapping of coefficients nreps times
    boot.coef <-  foreach(r=1:nreps, .combine='rbind') %dopar%  {

      ## extract new datasets from original dataset using the bootstrapped line numbers:
      boot.db <- db.modif [boot.ind$index[boot.ind$index$nrep == r, "ind" ], c("SpeciesCode", "abun",variable,'PlotName', covar)]

      # print(paste(r, ":","(",Sys.time(),")"))
      dat <- boot.db[as.character( boot.db$SpeciesCode)==sp, ] # select occurrences of the species
      names(dat) <- c("sp",'abun','var','PlotName', 'covar')[1:dim(dat)[2]]
      rm(boot.db)
      abun <- sort(as.numeric(as.character(na.omit(unique(dat$abun))))) ## list of abundance classes for species i

      ### Calculate GLM model if possible (sufficient abundance classes)
      if ("1" %in% abun & length(abun) >1 )
      {
        f <-  glm(dat$var ~ dat$covar + as.factor(dat$abun), family=poisson(log))

        # ## store coefficients for each bootstrapped sample k
        coef.string <- rep(NA,5)
        coef.covar <- summary(f)$coef [2, 1]
        for(j  in abun[-1]) coef.string[j-1] <- summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]

        if (!all(na.omit(coef.string) > -10)) coef.string <- rep(NA,5)
        print(paste(i, ":",r, ":", "(",Sys.time(),")"))
      }
       return(c(coef.covar, coef.string))
    }

    # assemble observed and bootstrapped coefs
    coefs <- rbind(c(covar.tab$coef.glm [i], obs.coef), boot.coef)
    covar.coefs <- coefs[,1]
    coefs <- coefs[,2:6]

    #### THIRD: statistics for species i
    
    # Bootstrapped coefs for the covariable
    covar.tab$mean[i] <- mean(covar.coefs, na.rm = T)
    covar.tab$sd[i] <- sd(covar.coefs, na.rm = T)
    covar.tab$lowCI[i] <-quantile(covar.coefs , probs = 0.025,na.rm =T)
    covar.tab$hiCI[i] <-quantile(covar.coefs , probs = 0.975,na.rm =T)

    # detect minimal critical value (simple start of negative trend)
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
    
  }

  stopCluster(cl)

  return(list(glms=glms,boot.mean = boot.mean, boot.sd=boot.sd, CIlow = bootCI.low, CIhi = bootCI.hi,
              covar.tab = covar.tab,
              crit.vals = crit.vals, spearman=spear, n.obs = n.obs,
              mean.values = mean.values,
              dif = dif,est= est, z=z,P= P))

}