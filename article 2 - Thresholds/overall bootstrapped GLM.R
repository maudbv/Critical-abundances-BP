### overall bootrapped GLM test : Elena's advice


###-----------Creating bootstrap samples of the dataset :
bootstrap.dataset<- function(db=databp[databp$PlotName %in% realgrasslands,], 
                           min.occur =5,  min.class = 1, R = 999) {
  
  ######### selecting species verifying conditions : 
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=1) 
                    &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
  db.modif <- db[which(db$SpeciesCode %in% a),] 
  
  # list of species =to be targeted in analysis :
  sp.names <- unique(db.modif$SpeciesCode)
  
  
  ######### Bootstrapping the dataset : resampling plots with replacement
    
    # Unique set of plots to be resampled 
    d = as.character(unique( db.modif$PlotName))
    
    # resampling the list of plot names with replacement *R* times
    boots <- sapply(1:R,FUN = function (k) sample(d, replace=TRUE))
    
    return(list(boots =boots,db.size=dim(db)[1],min.occur =min.occur,  min.class = min.class, R = R))
}
  

###-----------Calculating GLMs, critical values and impact size for bootstrapped samples

glm.overallboot<- function(db=databp[databp$PlotName %in% realgrasslands,], boot.output = boot.output,
                           variable='SR', min.occur =5,  min.class = 1, R = 999, call.boot=F) {
  
  # create bootstrap samples if necessary
  if (call.boot)  boot.output <- bootstrap.dataset(db=db, variable=variable, min.occur =min.occur,  min.class = min.class, R = R)
  
  ### Check that boot.output is the right one
  stopifnot(min.occur == boot.output$min.occur,
             min.class == boot.output$min.class,
             dim(db)[1] == boot.output$db.size)
  
  
  ######### selecting species verifying conditions : 
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=1) 
                    &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
  db.modif <- db[which(db$SpeciesCode %in% a),] 
  
  # list of species =to be targeted in analysis :
  sp.names <- unique(db.modif$SpeciesCode)
  
  
 #### Calculating GLM for observed and bootstrap datasets for each species

# Use the bootstrap samples
 boots<- boot.output$boots
 
# initiate result dataframes :
crit.vals <-data.frame(matrix(NA, nrow=length(sp.names), ncol= R +1,
                             dimnames=list(sp.names,1:(R+1))))
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

impact.size <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 8,
                             dimnames=list(sp.names, c("th","prevalence", "nb.plot.impact","prop.plot.impact",
                                                       "mean.dif","wtd.mean.dif", "th.dif", "max.dif"))))

impact.spread <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 4,
                                  dimnames=list(sp.names, c("th","prevalence", "nb.plot.impact","prop.plot.impact"))))


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
  
  coefs <- data.frame(matrix(NA, nrow=R+1, ncol= 5,
                            dimnames=list(1:(R+1) ,c("c2", "c3", "c4", "c5", "c6"))))
  
  for(j  in abun[-1]) coefs[1, j-1] = summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]

    #store coefs
    #   for(j  in abun[-1]) boot.coefs[[j-1]][1, k] = summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]
    #   

## SECOND : recalculate for each of the *R* bootstrapped datasets:
for ( r in 1:R) {
  
## extract new datasets from original dataset using the bootstrapped line numbers:
indices <-unlist(sapply(1:length(boots[,r]), function (k)  which(db.modif$PlotName == boots[k,r])))
boot.db <- db.modif [indices, c("SpeciesCode", "abun",variable,'PlotName')]  

print(paste(r, ":","(",Sys.time(),")"))
dat <- boot.db[as.character( boot.db$SpeciesCode)==sp, ] # select occurrences of the species
names(dat) <- c("sp",'abun','var','PlotName')[1:dim(dat)[2]]
abun <- sort(as.numeric(as.character(na.omit(unique(dat$abun))))) ## list of abundance classes for species i

### Calculate GLM model if possible (sufficient abundance classes)
if ("1" %in% abun & length(abun) >1 ) {
f <-  glm(dat$var ~ as.factor(dat$abun), family=poisson(log))

# ## store coefficients for each bootstrapped sample k
for(j  in abun[-1]) coefs[r+1, j-1] = summary(f)$coef [grep(paste(j), rownames(summary(f)$coef )), 1]
}

}

#### THIRD: statistics for species i 

# detect minimal critical value (simple start of negative trend)
crit.vals[i,] <- sapply(1:(R+1), function(x) {
crit <- NA
co <-coefs [x,]
ab = which(!is.na(co))+1
neg <-   which(co<0) +1

if (length(neg)>=1) {
  y <- neg [sapply(neg, FUN= function(l) {
    c1 <- ( if ( l+1 <= max(abun)) all(((l+1):max(ab)) %in% neg) # all higher classes have negative diferences
            else c1 =F)
    return(c1)
  })]
if (length(y)>=1) crit <- min( y, na.rm=T)
}
return(crit)
})

# calculate CI for bootstrapped coefs
boot.mean [i,] <-apply(coefs, 2, mean,na.rm =T) 
boot.sd [i,] <-apply(coefs, 2, sd,na.rm =T) 

bootCI.low [i,] <-apply(coefs, 2, quantile, probs = 0.025,na.rm =T) 
bootCI.hi [i,] <-apply(coefs, 2, quantile, probs = 0.975,na.rm =T) 

th <- crit.vals[i,1] 
prevalence <- dim(sp.dat)[1]
n.plot.impact<- sum(sp.dat$abun >=  crit.vals[i,1], na.rm=T)
prop.plot.impact <- n.plot.impact/prevalence
SRo <- mean(db[db$SpeciesCode == sp,variable], na.rm=T)
SRo <- mean(db[db$SpeciesCode == sp,variable], na.rm=T)

impact.spread[i,] <- c(th, prevalence, n.plot.impact, prop.plot.impact)
}

return(list(glms=glms, impact.spread = impact.spread,
            boot.mean = boot.mean, boot.sd=boot.sd, CIlow = bootCI.low, CIhi = bootCI.hi,
            crit.vals = crit.vals, spearman=spear, n.obs = n.obs,
            mean.values = mean.values,
            dif = dif,est= est, z=z,P= P))

}


### Calculate impact size and other indices 
## impact size for each species
impact.size <- function(obj = glmSR.overall){
  
  out = NULL
  
  for (i in 1:dim(obj$glms)[1]) {
    
  sp        <- rownames(obj$glms)[i] 
  th        <- obj$impact.spread$th[i]
  prevalence <- obj$impact.spread$prevalence[i]
  prop.plot.impact <- obj$impact.spread$prop.plot.impact[i]
  n.plot.impact <- obj$impact.spread$nb.plot.impact[i]
  
  dif = obj$dif[i,]
  n.obs = obj$dif[i,]
  SRo = obj$mean.values[i,"C1"]
  
  mean.dif <- wtd.mean.dif<-max.dif<- th.dif <- NA
  prop.mean.dif<- prop.wtd.mean.dif <-prop.max.dif <- impact.index <- NA
  
  #impact size for species with a threshold of impact :
  if (!is.na(th))  { 
    
    # mean impact
    mean.dif <- -mean(as.numeric(dif[c(th:6)-1]), na.rm=T)
    
    # frequency weighted mean impact
    d <- as.numeric(dif[c(th:6)-1])
    nb <- as.numeric(n.obs[c(th:6)-1])
    wtd.mean.dif <- -sum(d*nb, na.rm=T)/sum(nb, na.rm=T)
    
    # Max impact
    #! maximum dif is in fact the minimum because negative values
    max.dif <- - min(as.numeric(dif[c(th:6)-1]), na.rm=T) 
    
    # Threshold diference :
    th.dif = -as.numeric(dif[th-1])
    
    #########calculating proportional indices

    prop.mean.dif <- mean.dif /SRo
    prop.wtd.mean.dif <- wtd.mean.dif /SRo
    prop.max.dif <- max.dif /SRo
    
    ### overall impact index of the species
    impact.index <-prop.wtd.mean.dif * prop.plot.impact
  }
   out<- rbind(out, c(sp, th, prevalence, n.plot.impact, prop.plot.impact, 
                      mean.dif,  wtd.mean.dif,th.dif,max.dif,
                      SRo, prop.mean.dif, prop.wtd.mean.dif, prop.max.dif, impact.index))
  }
   out <- as.data.frame(out, as.is = T,stringsAsFactors = F)
   names(out) <- read.table(text = "sp,th,prevalence,n.plot.impact,prop.plot.impact,mean.dif,wtd.mean.dif,th.dif,max.dif,SRo,prop.mean.dif,prop.wtd.mean,prop.max.dif,impact.index",
                          stringsAsFactors = F, sep = ",")
  rownames(out) <- out$sp
return(out)
  }


