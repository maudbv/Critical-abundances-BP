
### Function recalculating the critical values using the bootstrapped CI
correct.impact.spread <- function(M = glmSRnat.overall, db=databp[databp$PlotName %in% unimprovedgrasslands,], variable = "SRnat", threshold.type = "custom") {

a <- row.names(M$dif)
db.modif <- db[which(db$SpeciesCode %in% a),]

# list of species =to be targeted in analysis :
sp.names <- a

output<- data.frame(matrix(NA, nrow= dim(M$dif), ncol =12))
names(output)<- c("th","th.CI.weak","th.CI","th.CI.custom", "pth" , "pth.CI.weak","pth.CI","pth.CI.custom","prevalence", "n.plot.impact","n.plot.dominant", "prop.plot.impact")
rownames(output) <- row.names(M$dif)

for (i in 1:length(sp.names)) {

  sp <- sp.names[i]  # select species name
  print(paste(i, ":", sp, "(",Sys.time(),")"))

  ### get GLM for observed dataset :
  sp.dat <- db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",variable,'PlotName')] # select occurrences of the species
  names(sp.dat) <- c('abun','var','PlotName')[1:dim(sp.dat)[2]]

# detect minimal critical value using bootstrap CI
  testCI = sign(as.numeric(M$CIlow [i,])) * sign(as.numeric(M$CIhi [i,]))
  co <-M$est [i,]
  ab = which(!is.na(co))+1
  ab.freq <- which(!is.na(co)& M$n.obs[i,2:6] >= min.occur) +1  ## with enough observations
  
  neg <-   which(co<0) +1

  th<- NA
  if (length(neg)>=1) {
    y <- neg [sapply(neg, FUN= function(l) {
      c1 <- ( if ( l <= max(ab)) all(((l):max(ab)) %in% neg) # all higher classes have negative diferences
              else c1 =F)
      return(c1)
    })]
    if (length(y)>=1) th <- min( y, na.rm=T)
  }

  # select coefficients whose CI are above or below zero :

  sig <-  which(testCI ==1 & M$n.obs[i,2:6] >= min.occur) + 1
  sig= sig [sig %in% neg]

# calculate critical abundance with all robust/significant negative coefs above
th.CI.weak<- NA
if (length(sig)>=1) {  # if there is at least  1 significant negative coef
    y <- sig [sapply(sig, FUN= function(l) {
      c1 <- ( if ( l <= max(ab)) all(((l):max(ab)) %in% neg) # all higher classes have negative diferences
              else c1 =F)
      return(c1)
    })]
    if (length(y)>=1) th.CI.weak <- min( y, na.rm=T)
}

th.CI<- NA
if (length(sig)>=1) {   # if there is at least  1 significant negative coef
  y <- sig [sapply(sig, FUN= function(l) { # for all significant coefs
    c1 <- ( if ( l <= max(ab)) (all(((l):max(ab.freq)) %in% sig) & all(((l):max(ab)) %in% neg)) # all higher classes have robust negative diferences
            else c1 =F)
    return(c1)
  })]
  if (length(y)>=1) th.CI <- min( y, na.rm=T)
}

th.CI.custom<- NA
if (length(sig)>=1) {   # if there is at least  1 significant negative coef
  y <- sig [sapply(sig, FUN= function(l) { # for all significant coefs
    c1 <- ( 
      if (l <= max(ab)) 
      {
         r1 = (all((l:max(ab)) %in% neg))  # all higher classes have negative coefficients, regardless of nobs.
        #r1 = (all((l:max(ab))[which(l:max(ab) %in% ab.freq)] %in% neg))  
        # all higher classes with sufficient nobs have negative coefficients
        
        r2 = all(((l):max(ab.freq)) %in% sig)
        # all higher classes have significantly significant (95% bootstrap CI) negative coefficients
        
        r3 = mean(as.numeric(M$est[i,l:max(ab)-1]), na.rm = T) <= M$est[i,l-1]
      # on average, native richness in higher classes are lower or equal than at the critical abundance
        
    c1 = r1 & (r2 | r3) # all negative coefficients AND significant OR on at least on average lower than the critical abundance
      }
            else 
              {
                c1 = FALSE
                }
      )
    return(c1)
  })]
  if (length(y)>=1) th.CI.custom <- min( y, na.rm=T)
}


### positive effects : detect minimal critical value using bootstrap CI
pos <- which(co>0) +1
# select coefficients whose CI are above or below zero :
sig <-  which(testCI ==1 & M$n.obs[i,2:6] >= min.occur) + 1
sig= sig [sig %in% pos]

pth<- NA
if (length(pos)>=1) {
  y <- pos [sapply(pos, FUN= function(l) {
    c1 <- ( if ( l <= max(ab)) all(((l):max(ab)) %in% pos) # all higher classes have negative diferences
            else c1 =F)
    return(c1)
  })]
  if (length(y)>=1) pth <- min( y, na.rm=T)
}

pth.CI.weak<- NA
if (length(sig)>=1) {  # if there is at least  1 significant negative coef
  y <- sig [sapply(sig, FUN= function(l) {
    c1 <- ( if ( l <= max(ab)) all(((l):max(ab)) %in% pos) # all higher classes have negative diferences
            else c1 =F)
    return(c1)
  })]
  if (length(y)>=1) pth.CI.weak <- min( y, na.rm=T)
}


pth.CI<- NA
if (length(sig)>=1) {
  y <- sig [sapply(sig, FUN= function(l) {
    c1 <- ( if ( l <= max(ab)) all(((l):max(ab.freq)) %in% sig) # all higher classes have negative diferences
            else c1 =F)
    return(c1)
  })]
  if (length(y)>=1) pth.CI <- min( y, na.rm=T)
}

pth.CI.custom<- NA
if (length(sig)>=1) {   # if there is at least  1 significant negative coef
  y <- sig [sapply(sig, FUN= function(l) { # for all significant coefs
    c1 <- ( 
      if (l <= max(ab)) 
      {
        r1 = (all((l:max(ab)) %in% pos))  # all higher classes have negative coefficients, regardless of nobs.
        #r1 = (all((l:max(ab))[which(l:max(ab) %in% ab.freq)] %in% neg))  
        # all higher classes with sufficient nobs have negative coefficients
        
        r2 = all(((l):max(ab.freq)) %in% sig)
        # all higher classes have significantly significant (95% bootstrap CI) negative coefficients
        
        r3 = mean(as.numeric(M$est[i,l:max(ab)-1]), na.rm = T) <= M$est[i,l-1]
        # on average, native richness in higher classes are lower or equal than at the critical abundance
        
        c1 = r1 & (r2 | r3) # all negative coefficients AND significant OR on at least on average lower than the critical abundance
      }
      else 
      {
        c1 = FALSE
      }
    )
    return(c1)
  })]
  if (length(y)>=1) pth.CI.custom <- min( y, na.rm=T)
}

## choosing type of threshold :

if (threshold.type == "weak") {
  th.CI <- th.CI.weak
  pth.CI <-pth.CI.weak
}


if (threshold.type == "custom") {
  th.CI <- th.CI.custom
  pth.CI <-pth.CI.custom
}

# Calculate impact spread

prevalence <- dim(sp.dat)[1]
nb.plot.impact<- sum( (sp.dat$abun >=  th.CI | sp.dat$abun >=  pth.CI ) , na.rm=T)
nb.plot.domin<- sum(sp.dat$abun >=  6, na.rm=T)
prop.plot.impact <- nb.plot.impact/prevalence

output [i, ] <- c(th, th.CI.weak, th.CI,th.CI.custom, pth, pth.CI.weak, pth.CI, pth.CI.custom, 
                  prevalence,nb.plot.impact, nb.plot.domin, prop.plot.impact)

rm(sig, pos, neg, th, th.CI.weak, th.CI,th.CI.custom, pth, pth.CI.weak, pth.CI, pth.CI.custom,
   prevalence,nb.plot.impact, nb.plot.domin, prop.plot.impact)

}


M$impact.spread <- output

return(M)
}
