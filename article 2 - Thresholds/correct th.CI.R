
### Function recalculating the critical values using the bootstrapped CI

correct.th.CI <- function(M = glmSR.overall, variable= "SR") {


a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=min.class) 
                  &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
db.modif <- db[which(db$SpeciesCode %in% a),] 

# list of species =to be targeted in analysis :
sp.names <- unique(db.modif$SpeciesCode)


for (i in 1:length(sp.names)) {
  
  sp <- sp.names[i]  # select species name
  print(paste(i, ":", sp, "(",Sys.time(),")"))
  
  ### FIRST : calculate GLM for observed dataset :
  sp.dat <- db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",variable,'PlotName')] # select occurrences of the species
  names(sp.dat) <- c('abun','var','PlotName')[1:dim(sp.dat)[2]]
  
# detect minimal critical value using bootstrap CI
  testCI = sign(as.numeric(M$CIlow [i,])) * sign(as.numeric(M$CIhi [i,]))
  co <-M$est [i,]
  ab = which(!is.na(co))+1
  neg <-   which(co<0) +1
  
  # select coefficients whose CI are above or below zero :
  sig <-  which(testCI ==1 & M$n.obs[i,2:6] >= min.occur) + 1
  sig= sig [sig %in% neg]
  
  th.CI<- NA
  if (length(sig)>=1) {
    y <- sig [sapply(sig, FUN= function(l) {
      c1 <- ( if ( l <= max(ab)) all(((l):max(ab)) %in% neg) # all higher classes have negative diferences
              else c1 =F)
      return(c1)
    })]
    if (length(y)>=1) th.CI <- min( y, na.rm=T)
}

M$impact.spread$th.CI [i] <- th.CI
M$impact.spread$nb.plot.impact [i] <- sum(sp.dat$abun >=  th.CI, na.rm=T)
M$impact.spread$prop.plot.impact [i] <- M$impact.spread$nb.plot.impact [i]/M$impact.spread$prevalence [i]
M$impact.spread$SRo [i] <- mean(db[db$SpeciesCode == sp,variable], na.rm=T)
}
return(M)
}
