### FUNCTIONS for GLM analyses of SR ~ abund classes

# fitting GLMs
#  = testing diference in a community metric (e.g. SR) between low abundances
# and increasing abundance classes using Generalized linear models 

glm.test <- function(db=databp[databp$PlotName %in% realgrasslands,], variable='SR', covar = NULL,
                     min.occur =5,  min.class = 1, alpha=0.05,CI=0.95, bootstrap = T, R = 999, drastic = F) {
  
  if (bootstrap) require(boot)
  
  
  ## sp occurring>=min.occur in the first abundance class
  ### and at least one abundance class of abundance 
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=1) 
           &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
  
  
  # select only species in a and b groups :
  db.modif <- db[which(db$SpeciesCode %in% a),] 
  
  # list of species name to be included in analysis :
  sp.names <- unique(db.modif$SpeciesCode)
  
  # initiate results
  init <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 5,
                             dimnames=list(sp.names,c("c2", "c3", "c4", "c5", "c6"))))
  dif <-   est<- P <- z <- init
  n.obs <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 6,
                            dimnames=list(sp.names,c("c1","c2", "c3", "c4", "c5", "c6"))))
  
  glms <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 3,
                           dimnames=list(sp.names, c( "df", "resid.dev","dev.ratio"))))
  
  spear <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 2,
                            dimnames=list(sp.names, c("rho" ,"p.val"))))
  
  thresh <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 9,
                             dimnames=list(sp.names, c("min.sig" ,"th",
                                                       "prevalence", "nb.plot.impact","prop.plot.impact",
                                                       "mean.dif","wtd.mean.dif", "th.dif", "max.dif"))))
  
  
  if (bootstrap)  {
    boot.table <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 6*6,
                                            dimnames=list(sp.names, 
                                            as.vector(sapply(1:6, FUN= function(k){
                                            paste("c", k, "_", c("t0","q2.5%","q97.5%","Pnegative","bias","se"), sep="")
                                             })))))
    boot.thresh  <-   thresh
  }
  
  # looping on species
  for (i in 1:length(sp.names) ) {
      
    sp <- sp.names[i]  # select species name
    
    print(paste(i, ":", sp, "(",Sys.time(),")"))
    
    sp.dat <- db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",variable,'PlotName', covar)] # select occurrences of the species
    names(sp.dat) <- c('abun','var','PlotName', "covar")[1:dim(sp.dat)[2]]
    
    abun <- sort(as.numeric(as.character(na.omit(unique(sp.dat$abun))))) ## list of abundance classes for species i
    n.obs[i,] <- sapply(1:6, FUN=function(j) length(sp.dat[which(sp.dat$abun==j),"var"]))
    
    
    #calculate diference in class mean SR
    dif[i,] <-  sapply(2:6, FUN=function(j) mean(sp.dat[which(sp.dat$abun==j),"var"], na.rm=T) - mean (sp.dat[which(sp.dat$abun==min(abun)),"var"], na.rm=T))
    
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
    
    ###### Bootstrapping   ###################################################
    boot.results <-  (if (bootstrap) {
    (fboot <- function(data = sp.dat) {   
      
    f.glm <- (if (!is.null(covar)) {
     function(d=sp.dat, w) { 
      f <-  glm(var ~ as.factor(abun) + covar ,data = d[w,], family=poisson(log))
      return(coefficients(f))
      }
    }
     else{
       function(d=sp.dat, w) { 
         f <-  glm(var ~ as.factor(abun),data = d[w,], family=poisson(log))
         return(coefficients(f))
       }
     }
     )
    
 boot.out <- boot(data, f.glm, R = R, strata = data$abun)   # stratifying by abundance level so that there is always at least one sample per abundance level
bootCI <- as.data.frame(t(sapply(1:(length(abun)+length(covar)), FUN =function(k) {
      bci <- boot.ci(boot.out, index = k,type=c("bca"), conf = CI)
      return(bci$bca)})))[,4:5]
    names(bootCI) = c("2.5%","97.5%")
    Pnegative <- apply(rbind(boot.out$t,boot.out$t0), 2, FUN= function(k) sum(k > 0)/(R+1))
    bias <-  colMeans(boot.out$t) - boot.out$t0
    boot.se <- sqrt(rowSums(apply(boot.out$t, 1, FUN= function(k) (k - colMeans(boot.out$t))^2)) / (R-1) )
    boot.results <- cbind(t0 = boot.out$t0, bootCI,  Pnegative= Pnegative, bias = bias, se = boot.se )
    
    

    return (boot.results)
      })()
    }
    else  NULL
    )
    
    boot.table [i, rep(1:6, length(abun)) + sort(rep((abun), 6)-1)*6 ] <- t(stack(as.data.frame(t(boot.results))))[1,]
   

#### Detect thresholds
    
    # which classes show negative dif
    neg <-   which(summary(f)$coef [, 3]<0) 
    
    ## which classes show signif negative dif
    sig  <-  which((summary(f)$coef [, 4] <= alpha) & ( summary(f)$coef [, 3] <0) & (n.obs[i,] >= min.occur)) 
    
    if (bootstrap) boot.sig  <-  which((sign(boot.results [1:length(abun),3]) == sign(boot.results [1:length(abun),2])) 
                                       & ( summary(f)$coef [1:length(abun), 3] <0)
                                       & (n.obs[i,] >= min.occur)) 


    ## alternative : 
    if (bootstrap & drastic) boot.sig  <-  which((sign(boot.results [1:length(abun),3]) == sign(boot.results [1:length(abun),2])) 
                                       & ( summary(f)$coef [1:length(abun), 3] <0
                                       & (summary(f)$coef [, 4] <= alpha))
                                       & (n.obs[i,] >= min.occur))
                                    


    #### Threshold detection
    thresh.detect <- function(signif, negative) {
    
    min.sig <- th <-  prevalence <- n.plot.impact <- mean.dif <- wtd.mean.dif <- th.dif <- max.dif <- NA
    
    ### identifying a threshold 
    ###  => lowest abundance where significant negative impact, followed by only negative trends (significant or not).
    if (length(signif)!=0) {
      min.sig <-  min(signif)    # init threshold value with lowest class with signif negative diference
      test <-  F
      j <- 0
      
      # search for thresholds in case min.sig is not followed by negative diferences
      th.sig = sapply(signif, FUN= function(k) {
        c1 <- all(((k):max(abun)) %in% negative) # all higher classes have negative diferences
        c2 <- (length(sp.dat[sp.dat$abun == k, "var"]) >= min.occur) # at least min.occur observation in the threshold class
        return(c1 & c2)
      })
      
      if (! all(th.sig == F)) th = min(signif[th.sig == T])
    }     
    
     
    ## additional stats per species :
    prevalence <- dim(sp.dat)[1]
    n.plot.impact <- sum(sp.dat$abun >= th, na.rm=T)
    prop.plot.impact <- n.plot.impact/prevalence
    
    #impact size for species with a threshold of impact :
    if (!is.na(th))  { 
    
        # mean impact
        mean.dif <- -mean(as.numeric(dif[i,c(th:6)-1]), na.rm=T)
        
        # frequency weighted mean impact
        d <- as.numeric(dif[i, c(th:6)-1])
        nb <- as.numeric(n.obs[i, c(th:6)-1])
        wtd.mean.dif <- -sum(d*nb, na.rm=T)/sum(nb, na.rm=T)
        
        # Max impact
        #! maximum dif is in fact the minimum because negative values
        max.dif <- - min(as.numeric(dif[i, c(th:6)-1]), na.rm=T) 
  
        # Threshold diference :
        th.dif = -as.numeric(dif[i,th-1])
    }
    
  return( c(min.sig, th, prevalence, n.plot.impact, prop.plot.impact, mean.dif,  wtd.mean.dif,th.dif,max.dif))
    }

thresh[i,] <- thresh.detect(signif = sig, negative = neg)
boot.thresh[i,] <- thresh.detect(signif = boot.sig, negative = neg)
  }
  
if (bootstrap) return(list(glms=glms, boot= boot.table, thresh = thresh, boot.thresh= boot.thresh, spearman=spear, n.obs = n.obs,  dif = dif,est= est, z=z,P= P))
if (!bootstrap) return(list(glms=glms, thresh = thresh,  spearman=spear, n.obs = n.obs,  dif = dif,est= est, z=z,P= P))

}

# fitting GLMs
#  and bootstrapping with no strata : TO BE DONE

glm.test.nostrata <- function(db=databp[databp$PlotName %in% realgrasslands,], variable='SR', covar = NULL,
                     min.occur =5,  min.class = 1, alpha=0.05,CI=0.95, bootstrap = T, R = 999, drastic = F) {
  
  if (bootstrap) require(boot)
  
  
  ## sp occurring>=min.occur in the first abundance class
  ### and at least one abundance class of abundance 
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=1) 
                    &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
  
  
  # select only species in a and b groups :
  db.modif <- db[which(db$SpeciesCode %in% a),] 
  
  # list of species name to be included in analysis :
  sp.names <- unique(db.modif$SpeciesCode)
  
  # initiate results
  init <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 5,
                             dimnames=list(sp.names,c("c2", "c3", "c4", "c5", "c6"))))
  dif <-   est<- P <- z <- init
  n.obs <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 6,
                            dimnames=list(sp.names,c("c1","c2", "c3", "c4", "c5", "c6"))))
  
  glms <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 3,
                           dimnames=list(sp.names, c( "df", "resid.dev","dev.ratio"))))
  
  spear <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 2,
                            dimnames=list(sp.names, c("rho" ,"p.val"))))
  
  thresh <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 9,
                               dimnames=list(sp.names, c("min.sig" ,"th",
                                                         "prevalence", "nb.plot.impact","prop.plot.impact",
                                                         "mean.dif","wtd.mean.dif", "th.dif", "max.dif"))))
  
  
  if (bootstrap)  {
    boot.table <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 6*6,
                                   dimnames=list(sp.names, 
                                                 as.vector(sapply(1:6, FUN= function(k){
                                                   paste("c", k, "_", c("t0","q2.5%","q97.5%","Pnegative","bias","se"), sep="")
                                                 })))))
    boot.thresh  <-   thresh
  }
  
  # looping on species
  for (i in 1:length(sp.names) ) {
    
    sp <- sp.names[i]  # select species name
    
    print(paste(i, ":", sp, "(",Sys.time(),")"))
    
    sp.dat <- db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",variable,'PlotName', covar)] # select occurrences of the species
    names(sp.dat) <- c('abun','var','PlotName', "covar")[1:dim(sp.dat)[2]]
    
    abun <- sort(as.numeric(as.character(na.omit(unique(sp.dat$abun))))) ## list of abundance classes for species i
    n.obs[i,] <- sapply(1:6, FUN=function(j) length(sp.dat[which(sp.dat$abun==j),"var"]))
    
    
    #calculate diference in class mean SR
    dif[i,] <-  sapply(2:6, FUN=function(j) mean(sp.dat[which(sp.dat$abun==j),"var"], na.rm=T) - mean (sp.dat[which(sp.dat$abun==min(abun)),"var"], na.rm=T))
    
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
    
    ###### Bootstrapping   ###################################################
    boot.results <-  (if (bootstrap) {
      (fboot <- function(data = sp.dat) {   
        
        f.glm <- function(d=sp.dat, w) { 
            f <-  glm(var ~ as.factor(abun),data = d[w,], family=poisson(log))
            return(coefficients(f))
          }

        boot.out <- boot(data, f.glm, R = R)   # no stratifying : sample size changes
        bootCI <- as.data.frame(t(sapply(1:(length(abun)+length(covar)), FUN =function(k) {
          bci <- boot.ci(boot.out, index = k,type=c("bca"), conf = CI)
          return(bci$bca)})))[,4:5]
        names(bootCI) = c("2.5%","97.5%")
        Pnegative <- apply(rbind(boot.out$t,boot.out$t0), 2, FUN= function(k) sum(k > 0)/(R+1))
        bias <-  colMeans(boot.out$t) - boot.out$t0
        boot.se <- sqrt(rowSums(apply(boot.out$t, 1, FUN= function(k) (k - colMeans(boot.out$t))^2)) / (R-1) )
        boot.results <- cbind(t0 = boot.out$t0, bootCI,  Pnegative= Pnegative, bias = bias, se = boot.se )
        
        
        
        return (boot.results)
      })()
    }
    else  NULL
    )
    
    boot.table [i, rep(1:6, length(abun)) + sort(rep((abun), 6)-1)*6 ] <- t(stack(as.data.frame(t(boot.results))))[1,]
    
    
    #### Detect thresholds
    
    # which classes show negative dif
    neg <-   which(summary(f)$coef [, 3]<0) 
    
    ## which classes show signif negative dif
    sig  <-  which((summary(f)$coef [, 4] <= alpha) & ( summary(f)$coef [, 3] <0) & (n.obs[i,] >= min.occur)) 
    
    if (bootstrap) boot.sig  <-  which((sign(boot.results [1:length(abun),3]) == sign(boot.results [1:length(abun),2])) 
                                       & ( summary(f)$coef [1:length(abun), 3] <0)
                                       & (n.obs[i,] >= min.occur)) 
    
    
    ## alternative : 
    if (bootstrap & drastic) boot.sig  <-  which((sign(boot.results [1:length(abun),3]) == sign(boot.results [1:length(abun),2])) 
                                                 & ( summary(f)$coef [1:length(abun), 3] <0
                                                     & (summary(f)$coef [, 4] <= alpha))
                                                 & (n.obs[i,] >= min.occur))
    
    
    
    #### Threshold detection
    thresh.detect <- function(signif, negative) {
      
      min.sig <- th <-  prevalence <- n.plot.impact <- mean.dif <- wtd.mean.dif <- th.dif <- max.dif <- NA
      
      ### identifying a threshold 
      ###  => lowest abundance where significant negative impact, followed by only negative trends (significant or not).
      if (length(signif)!=0) {
        min.sig <-  min(signif)    # init threshold value with lowest class with signif negative diference
        test <-  F
        j <- 0
        
        # search for thresholds in case min.sig is not followed by negative diferences
        th.sig = sapply(signif, FUN= function(k) {
          c1 <- all(((k):max(abun)) %in% negative) # all higher classes have negative diferences
          c2 <- (length(sp.dat[sp.dat$abun == k, "var"]) >= min.occur) # at least min.occur observation in the threshold class
          return(c1 & c2)
        })
        
        if (! all(th.sig == F)) th = min(signif[th.sig == T])
      }     
      
      
      ## additional stats per species :
      prevalence <- dim(sp.dat)[1]
      n.plot.impact <- sum(sp.dat$abun >= th, na.rm=T)
      prop.plot.impact <- n.plot.impact/prevalence
      
      #impact size for species with a threshold of impact :
      if (!is.na(th))  { 
        
        # mean impact
        mean.dif <- -mean(as.numeric(dif[i,c(th:6)-1]), na.rm=T)
        
        # frequency weighted mean impact
        d <- as.numeric(dif[i, c(th:6)-1])
        nb <- as.numeric(n.obs[i, c(th:6)-1])
        wtd.mean.dif <- -sum(d*nb, na.rm=T)/sum(nb, na.rm=T)
        
        # Max impact
        #! maximum dif is in fact the minimum because negative values
        max.dif <- - min(as.numeric(dif[i, c(th:6)-1]), na.rm=T) 
        
        # Threshold diference :
        th.dif = -as.numeric(dif[i,th-1])
      }
      
      return( c(min.sig, th, prevalence, n.plot.impact, prop.plot.impact, mean.dif,  wtd.mean.dif,th.dif,max.dif))
    }
    
    thresh[i,] <- thresh.detect(signif = sig, negative = neg)
    boot.thresh[i,] <- thresh.detect(signif = boot.sig, negative = neg)
  }
  
  if (bootstrap) return(list(glms=glms, boot= boot.table, thresh = thresh, boot.thresh= boot.thresh, spearman=spear, n.obs = n.obs,  dif = dif,est= est, z=z,P= P))
  if (!bootstrap) return(list(glms=glms, thresh = thresh,  spearman=spear, n.obs = n.obs,  dif = dif,est= est, z=z,P= P))
  
}

#########calculating proportional indices
calc.prop.impact = function(variable="SR" , M = glmSR, db = db , index = "mean.dif") {
  X <-M$boot.thresh
  # X <-M$boot.thresh[!is.na(M$boot.thresh$th),]
  Z<- db[db$SpeciesCode %in% rownames(X),]
  tab <- tapply(Z[,variable], list(Z$SpeciesCode , Z$abun), FUN=mean, na.rm=T)
  prop.index <- X[,index] /tab[,1]
  return(prop.index)
}

      # adding prop impact size to the table
  add.prop <- function(N = glmSR, variable ="SR", data=db) { 
        modif <-N$boot.thresh
        modif$prop.mean.dif <- calc.prop.impact(variable="SR" ,M=N, db = data, index = "mean.dif")
        modif$prop.wtd.mean.dif <- calc.prop.impact(variable="SR" , M =N, db = data, index = "wtd.mean.dif")
        modif$prop.max.dif <- calc.prop.impact(variable="SR" , M =N, db = data, index = "max.dif")
        
        modif$index <-sqrt(modif$prop.wtd.mean.dif) * sqrt(modif$prop.plot.impact)
        
        return(modif)
      }

# summarizing results on thresholds per group of species
summary.glmtest <- function(M = glmSR.grass,data=species, group="ALIEN", type =c("glm", "boot")) {
  
  ### select onlys species in the model
  data=data[rownames(M$dif),]
  
   if (type == "boot") M$thresh <- M$boot.thresh
  
  ## divide species according to grouping factor (ie Alien vs Native)
  if (is.null(group)) {
    sub<-M
    G<-1
  }
  
  if (!is.null(group)) { 
    sub<-list()
    i<-1
    G<-length(unique(data[,group]))
    for (g in sort(unique(data[,group]))) {
      sub[[i]]<-lapply(M, FUN= function(x) x[rownames(x) %in% rownames(data)[data[,group]==g],] )
      i<-i+1
    }
    names<-paste(group, ":",sort(unique(data[,group])), sep="")
  }                                 
  
  out<-NULL   
  sum.table <- NULL 
  
  for (j in 1:G) {
    S<-sub[[j]]
    
    # number of potentieally signif species occurrence per class
    n.sp<-sapply(names(S$dif), FUN= function(x) {
      length(which(!is.na(S$dif[,x])))
    })
    
    # number of species which have a threshold per class
    n.target<-sapply(names(S$dif), FUN= function(x) {
      length(which(!is.na(S$dif[,x]) & !is.na(S$thresh$th)) )
    })
    
    
    # number of significant negative effects per class
    n.negative<-sapply(names(S$dif), FUN= function(x) {
      length(which(S$P[,x]<=0.05 & S$z[,x]<0))
    })
    
    
    # number of significant positive effects per class
    n.positive<-sapply(names(S$dif), FUN= function(x) {
      length(which(S$P[,x]<0.05 & S$z[,x]>=0))
    })
    
    # number of significant negative impavts per class for species with a threshld
    n.impact<-sapply(names(S$dif), FUN= function(x) {
      length(which(S$P[,x]<=0.05 & S$z[,x]<0 & !is.na(S$thresh$th)))
    })
    
    # proportion of signifi negative effects per class
    p.impact<-n.impact/ n.sp
    
    # proportion of signifi negative effects per class
    p.negative<-n.negative/ n.sp
    
    # number of times each class is the threshold
    freq.thr <-table(as.factor(S$thresh$th))[match(2:6,names(table(as.factor(S$thresh$th))))]
    names(freq.thr)<-names(n.sp)
    freq.thr[is.na(freq.thr)]<-0
    prop.thr <- freq.thr / n.sp
    if(sum(freq.thr, na.rm=T)!=0) perc.thr <- freq.thr / sum(freq.thr, na.rm=T)
    if(sum(freq.thr, na.rm=T)==0) perc.thr <- freq.thr 
    
    out<-rbind(out, data.frame(group=names[j],class = 2:6, nb.sp=n.sp, n.target = n.target, freq.impact=n.impact, freq.negative =n.negative,freq.positive =n.positive,
                              prop.impact=as.array(p.impact),freq.thr=freq.thr, prop.thr=prop.thr, perc.thr=perc.thr))
    
    #summary table
    n.sp <- length(rownames(S$thresh))
    freq.thr <- length(which(!is.na(S$thresh$th)))
       
    sum.table <-rbind(sum.table,  data.frame(group=names[j], nb.sp=n.sp,  freq.thr=freq.thr))

  }
  return(list(class.summary =out, overall.summary = sum.table))
}

