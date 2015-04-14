
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

