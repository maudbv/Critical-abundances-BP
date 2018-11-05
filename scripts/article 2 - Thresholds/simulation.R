#simulation of random species loss

# reference community matrix:
ref.comm <- comm[ unimprovedgrasslands, names(comm) %in% natives]
ref.comm <- ref.comm [rowSums(ref.comm, na.rm = T) != 0,colSums(ref.comm, na.rm = T) != 0]
dim(ref.comm)

S <- dim(ref.comm)[2] # number of native species in the gamma pool
N <- dim(ref.comm)[1] # number of plots
R <- 20 # mean native richness 


# Extract frequency distribution of native species across the landscape based on our data  in BP :
freq <- colSums(ref.comm>0, na.rm = T)
rel.freq  <- freq/N

## Frequency of each species per abundance class when present
freq.wp <- t(apply(ref.comm, 2, function(x) table(factor(x, levels = c(0,1,2,3,4,5,6)))))[,-1]
m.sp <-  apply(freq.wp, 2, mean)
s.sp <- apply(freq.wp, 2, sd)/2  ### /dim(y)[2]


## Abudance distribution within native communities in data:
sad <- t(apply(ref.comm, 1, function(x) table(factor(x, levels = c(0,1,2,3,4,5,6)))))
m.sad <-  apply(sad[,-1], 2, mean)
s.sad <- apply(sad[,-1], 2, sd)/2  ### /dim(y)[2]

## Check also alien distribtuions:
ali.comm <- comm[ unimprovedgrasslands, names(comm) %in% aliens]
ali.comm <- ali.com[rowSums(ali.comm, na.rm = T) != 0,colSums(ali.comm, na.rm = T) != 0]

sad.ali <- t(apply(ali.comm, 1, function(x) table(factor(x, levels = c(0,1,2,3,4,5,6)))))
m.sad.ali <-  apply(sad.ali[,-1], 2, mean)
s.sad.ali <- apply(sad.ali[,-1], 2, sd)/2  #/dim(y.ali)[1]


par(mfrow = c(1,2))
b <- barplot(m.sad, ylim = c(0,6))
segments(b, m.sad-s.sad, b, m.sad+s.sad)
b <- barplot(m.sad.ali, ylim = c(0,10))
segments(b, m.sad.ali-s.sad.ali, b, m.sad.ali + s.sad.ali)



## SETUP community matrix with only presence/absence: BOF ############
SRnat <- rpois(N, R) 

r.comm <- t(sapply(1:N, function(i) {
  x = rep(0,S) 
  x [sample( 1:S, SRnat[i], prob = as.numeric(freq))] <- 1
  return(x)
}))

##  SIMULATE INITIAL COMMUNITY FAILED  #######

# attribute natives richness for each of the N plots so that it is on average R species
SRnat <- rpois(N, R) 

## Simulate abundance distributions for natives:
sad.sim <- as.matrix(data.frame("1" = rpois(N, m.sad[1]),
                                "2" = rpois(N, m.sad[2]), 
                                "3" = rpois(N, m.sad[3]), 
                                "4" = rpois(N, m.sad[4]), 
                                "5" = rpois(N, m.sad[5]),
                                "6" = rpois(N, m.sad[6])
))

m.sad.sim <-  apply(sad.sim, 2, mean)
s.sad.sim <- apply(sad.sim, 2, sd)/2

## Compare simulated SAD to observed SAD :
par(mfrow = c(1,2))
b <- barplot(m.sad, ylim = c(0,5))
segments(b, m.sad-s.sad, b, m.sad+s.sad)

b <- barplot(m.sad.sim, ylim = c(0,5))
segments(b, m.sad.sim-s.sad.sim, b, m.sad.sim+s.sad.sim)

#Check visually random sets of distributions:
par(mfrow = c(1,2))
i = ceiling(runif(1,1,N))
barplot(sad[i,-1] ,  ylim = c(0,10),main = "observed")
barplot(sad.sim[i,], ylim = c(0,10), main ="simulated")

## Distribution of local native richness 
hist(rowSums(sad[,-1]))  # typical poisson : mean = sd
hist(rowSums(sad.sim))   # typical Gaussian : same mean as oberved, but lower sd

### PROBLEM with SIMULATION of SAD !! (not surprisingly...) ##


## Simulate distribution of abundance classes for each species:
freq.wp.sim <- as.matrix(data.frame("1" = rpois(S, m.sp[1]),
                                    "2" = rpois(S, m.sp[2]), 
                                    "3" = rpois(S, m.sp[3]), 
                                    "4" = rpois(S, m.sp[4]), 
                                    "5" = rpois(S, m.sp[5]),
                                    "6" = rpois(S, m.sp[6])
))

m.sp.sim <-  apply(freq.wp.sim, 2, mean)
se.sp.sim <- apply(freq.wp.sim, 2, sd)/2  ### /dim(y)[2]
rel.freq.wp.sim <- freq.wp.sim/rowSums(freq.wp.sim)

## Compare simulated SAD to observed SAD :
par(mfrow = c(1,2))
b <- barplot(m.sp, ylim = c(0,5))
segments(b, m.sp-s.sp, b, m.sp+s.sp)
b <- barplot(m.sp.sim, ylim = c(0,6))
segments(b, m.sp.sim-se.sp.sim, b, m.sp.sim+se.sp.sim)


## Distribution of species which are always rare at local abundance   
# == more than 80 % of the time either rare or occasional
par(mfrow = c(1,2))
rarity <- rowSums(freq.wp[,1:2]) * (rowSums(rel.freq.wp[,1:2])>0.8)
barplot(table(rarity))
rarity.sim <- rowSums(freq.wp.sim[,1:2]) * (rowSums(rel.freq.wp.sim[,1:2])>0.8)
barplot(table(rarity.sim))
## PROBLEM : 
# probably the frequency of species in each abudnance class is not a poisson distribution, 
# but more a conditional binomial ? 

rare.sp <- names(rarity)[which(rarity > 0)]

# ====> DROP IT




## SETUP  with real observed values  #######
r.comm <- ref.comm
SRnat <- S

## Calculate initial values

g.below = sum(colSums(r.comm[1:N/2,])>0)
g.above = sum(colSums(r.comm[(N/2 + 1):N,])>0)
Dg.init <- g.below-g.above

a.below = mean(rowSums(r.comm[1:N/2,]))
a.above = mean(rowSums(r.comm[(N/2 + 1):N,]))
Da.init <-  a.below - a.above

### SIMULATE LOSS of always the SAME species:  #####

nreps = 50
max.loss = 20 

# random species
simul.rand.same <- t( sapply(0:max.loss, function(l) {
  runs <- as.data.frame(t( sapply(1:nreps, function(rep) {
    sp <- sample(1:S, l) # random species, but always the same lost
    mod.comm <- r.comm
    mod.comm[(N/2 + 1):N, sp] <- 0   # always the same  species disappearing in all plots above
    
    g.below = sum(colSums(mod.comm[1:N/2,])>0)
    g.above = sum(colSums(mod.comm[(N/2 + 1):N,])>0)
    Dg <- g.below-g.above
    
    a.below = mean(rowSums(mod.comm[1:N/2,]))
    a.above = mean(rowSums(mod.comm[(N/2 + 1):N,]))
    Da <-  a.below-a.above
    
    return( c(Da= Da, Dg = Dg))
  })))
  print(l)
  return( c(loss= l, nreps = rep, Da = mean(runs$Da), Da.sd = sd (runs$Da), Dg = mean(runs$Dg), Dg.sd = sd(runs$Dg)) )
}))

plot(Dg ~ loss, simul.rand.same, ylim =c(0,30))
points(Da ~ loss, simul.rand.same, col = "red")
plot(Dg ~ Da, simul.rand.same)

# Only rare species, always the same
simul.rare.same <- t( sapply(1:30, function(l) {
  runs <- as.data.frame(t( sapply(1:nreps, function(rep) {
    sp <- sample(1:S, l,  prob = as.numeric(1 - (freq))) # rare species are more likely to disappear locally
    mod.comm <- r.comm
    mod.comm[(N/2 + 1):N, sp] <- 0   # always the same  species disappearing in all plots above
    g.below = sum(colSums(mod.comm[1:N/2,])>0)
    g.above = sum(colSums(mod.comm[(N/2 + 1):N,])>0)
    Dg <- g.below-g.above
    Dg
    
    a.below = mean(rowSums(mod.comm[1:N/2,]))
    a.above = mean(rowSums(mod.comm[(N/2 + 1):N,]))
    Da <-  a.below-a.above
    Da
    
    return( c(Da= Da, Dg = Dg))
  })))
  print(l)
  return( c(loss= l, nreps = rep, Da = mean(runs$Da), Da.sd = sd (runs$Da), Dg = mean(runs$Dg), Dg.sd = sd(runs$Dg)) )
}))

plot(Dg ~ loss, simul.rare.same, ylim =c(0,30))
points(Da ~ loss, simul.rare.same, col = "red")
plot(Dg ~ Da, simul.rare.same)


# Only freq species, always the same
simul.freq.same <- t( sapply(0:30, function(l) {
  runs <- as.data.frame(t( sapply(1:nreps, function(rep) {
    sp <- sample(1:S, l,  prob = as.numeric(freq)) # most freq species are more likely to disappear locally
    mod.comm <- r.comm
    mod.comm[(N/2 + 1):N, sp] <- 0   # always the same  species disappearing in all plots above
    
    g.below = sum(colSums(mod.comm[1:N/2,])>0)
    g.above = sum(colSums(mod.comm[(N/2 + 1):N,])>0)
    Dg <- g.below-g.above
    Dg
    
    a.below = mean(rowSums(mod.comm[1:N/2,]))
    a.above = mean(rowSums(mod.comm[(N/2 + 1):N,]))
    Da <-  a.below-a.above
    Da
    
    return( c(Da= Da, Dg = Dg))
  })))
  print(l)
  return( c(loss= l, nreps = rep, Da = mean(runs$Da), Da.sd = sd (runs$Da), Dg = mean(runs$Dg), Dg.sd = sd(runs$Dg)) )
}))

plot(Dg ~ loss, simul.freq.same, ylim =c(0,40))
points(Da ~ loss, simul.freq.same, col = "red")
plot(Dg ~ Da, simul.freq.same)


## Loss of DIFFERENT species in each plot  #######

# random species
simul.rand.diff <- t( sapply(1:20, function(l) {
  runs <- as.data.frame(t( sapply(1:nreps, function(rep) {
    
    mod.comm <- r.comm
    mod.comm[(N/2 + 1):N,] <-  t(apply(mod.comm[(N/2 + 1):N,], 1, function(x) {
      pres <- which(x >0)
      x[as.numeric(sample(as.character(pres), min(length(pres),rpois(1,l))))] <- 0
      return(x)
    }))
    
    g.below = sum(colSums(mod.comm[1:N/2,])>0)
    g.above = sum(colSums(mod.comm[(N/2 + 1):N,])>0)
    Dg <- g.below-g.above
    Dg
    
    a.below = mean(rowSums(mod.comm[1:N/2,]))
    a.above = mean(rowSums(mod.comm[(N/2 + 1):N,]))
    Da <-  a.below-a.above
    Da
    
    return( c(l = l , Da= Da, Dg = Dg))
  })))
  print(l)
  return( c(loss= l , nreps = rep, Da = mean(runs$Da), Da.sd = sd (runs$Da), Dg = mean(runs$Dg), Dg.sd = sd(runs$Dg)) )
}))

par( mfrow = c(1,3))
plot(Dg ~ loss, simul.rand.diff)
abline(Dg.init,1)
plot(Da ~ loss, simul.rand.diff, col = "red", ylim = c(0,20))
abline(Da.init,1)
plot(Dg ~ Da, simul.rand.diff)
abline(Dg.init,1)

# Locally rare species
simul.rare.diff <- t( sapply(1:20, function(l) {
  runs <- as.data.frame(t( sapply(1:nreps, function(rep) {
    
    mod.comm <- r.comm
    mod.comm[(N/2 + 1):N,] <-  t(apply(mod.comm[(N/2 + 1):N,], 1, function(x) {
      pres <- which(x %in% c(1,2))  ## select only among the locally rare or occasional species
      x[as.numeric(sample(as.character(pres), min(length(pres),rpois(1,l))))] <- 0
      return(x)
    }))
    
    g.below = sum(colSums(mod.comm[1:N/2,])>0)
    g.above = sum(colSums(mod.comm[(N/2 + 1):N,])>0)
    Dg <- g.below-g.above
    Dg
    
    a.below = mean(rowSums(mod.comm[1:N/2,]))
    a.above = mean(rowSums(mod.comm[(N/2 + 1):N,]))
    Da <-  a.below-a.above
    Da
    
    return( c(l = l , Da= Da, Dg = Dg))
  })))
  print(l)
  return( c(loss= l , nreps = rep, Da = mean(runs$Da), Da.sd = sd (runs$Da), Dg = mean(runs$Dg), Dg.sd = sd(runs$Dg)) )
}))


par( mfrow = c(1,3))
plot(Dg ~ loss, simul.rare.diff)
abline(Dg.init,1)
plot(Da ~ loss, simul.rare.diff, col = "red", ylim = c(0,20))
abline(Da.init,1)
plot(Dg ~ Da, simul.rare.diff)
abline(Dg.init,1)


# INFREQUENT  species
# version with all frequencies
simul.infr.diff <- t( sapply(1:20, function(l) {
  runs <- as.data.frame(t( sapply(1:nreps, function(rep) {

    mod.comm <- r.comm
    mod.comm[(N/2 + 1):N,] <-  t(apply(mod.comm[(N/2 + 1):N,], 1, function(x) {
      pres = which(x>0)
      x[as.numeric(sample(as.character(pres), min(length(pres),rpois(1,l)),
                          prob = as.numeric(1 - (rel.freq[pres]))
                          ))] <- 0
      return(x)
    }))
    
    g.below = sum(colSums(mod.comm[1:N/2,])>0)
    g.above = sum(colSums(mod.comm[(N/2 + 1):N,])>0)
    Dg <- g.below-g.above
    Dg
    
    a.below = mean(rowSums(mod.comm[1:N/2,]))
    a.above = mean(rowSums(mod.comm[(N/2 + 1):N,]))
    Da <-  a.below-a.above
    Da
    
    return( c(l = l , Da= Da, Dg = Dg))
  })))
  print(l)
  return( c(loss= l , nreps = rep, Da = mean(runs$Da), Da.sd = sd (runs$Da), Dg = mean(runs$Dg), Dg.sd = sd(runs$Dg)) )
}))


# version with only species < 10 plots
table(freq < 10 )

simul.infr.diff <- t( sapply(1:20, function(l) {
  runs <- as.data.frame(t( sapply(1:nreps, function(rep) {
    
    mod.comm <- r.comm
    mod.comm[(N/2 + 1):N,] <-  t(apply(mod.comm[(N/2 + 1):N,], 1, function(x) {
      pres = which(x>0) 
      pres <- pres[freq[pres] <10]   ## Keep only the infrequent (less than 10 plots)
      x[as.numeric(sample(as.character(pres), min(length(pres),rpois(1,l))))] <- 0
      return(x)
    }))
    
    g.below = sum(colSums(mod.comm[1:N/2,])>0)
    g.above = sum(colSums(mod.comm[(N/2 + 1):N,])>0)
    Dg <- g.below-g.above
    Dg
    
    a.below = mean(rowSums(mod.comm[1:N/2,]))
    a.above = mean(rowSums(mod.comm[(N/2 + 1):N,]))
    Da <-  a.below-a.above
    Da
    
    return( c(l = l , Da= Da, Dg = Dg))
  })))
  print(l)
  return( c(loss= l , nreps = rep, Da = mean(runs$Da), Da.sd = sd (runs$Da), Dg = mean(runs$Dg), Dg.sd = sd(runs$Dg)) )
}))

par( mfrow = c(1,3))
plot(Dg ~ loss, simul.infr.diff)
abline(Dg.init,1)
plot(Da ~ loss, simul.infr.diff, col = "red", ylim = c(0,20))
abline(Da.init,1)
plot(Dg ~ Da, simul.infr.diff)
abline(Dg.init,1)

#freq species
simul.freq.diff <- t( sapply(1:20, function(l) {
  runs <- as.data.frame(t( sapply(1:nreps, function(rep) {
    
    mod.comm <- r.comm
    mod.comm[(N/2 + 1):N,] <-  t(apply(mod.comm[(N/2 + 1):N,], 1, function(x) {
      pres = which(x>0)
      x[as.numeric(sample(as.character(pres), min(length(pres),rpois(1,l)),
                          prob = as.numeric(rel.freq[pres])
      ))] <- 0
      return(x)
    }))
    
    g.below = sum(colSums(mod.comm[1:N/2,])>0)
    g.above = sum(colSums(mod.comm[(N/2 + 1):N,])>0)
    Dg <- g.below-g.above
    Dg
    
    a.below = mean(rowSums(mod.comm[1:N/2,]))
    a.above = mean(rowSums(mod.comm[(N/2 + 1):N,]))
    Da <-  a.below-a.above
    Da
    
    return( c(l = l , Da= Da, Dg = Dg))
  })))
  print(l)
  return( c(loss= l , nreps = rep, Da = mean(runs$Da), Da.sd = sd (runs$Da), Dg = mean(runs$Dg), Dg.sd = sd(runs$Dg)) )
}))

plot(Dg ~ loss, simul.freq.diff, ylim =c(0,100))
points(Da ~ loss, simul.freq.diff, col = "red")
plot(Dg ~ Da, simul.freq.diff)


## Summary graphes  #####
plot(Dg ~ loss, simul.rare.diff, ylim =c(Dg.init,100), type = "l",col = "red")
lines(Dg ~ loss, simul.rand.diff,col = "blue")
lines(Dg ~ loss, simul.freq.diff, col = "forestgreen")
lines(Dg ~ loss, simul.infr.diff, col = "violet")

par(mfrow = c(1,1))
plot(Dg ~ Da, simul.rare.diff, ylim =c(Dg.init,90), xlim = c(Da.init,14), type = "l",col = "red")
lines(Dg ~ Da, simul.rand.diff,col = "black")
lines(Dg ~ Da, simul.rand.same,col = "black", lty = "dashed")

lines(Dg ~ Da, simul.freq.diff, col = "forestgreen")
lines(Dg ~ Da, simul.freq.same, col = "forestgreen", lty = "dashed")

lines (Dg ~ Da, simul.rare.same, col = "violet", lty = "dashed")
lines (Dg ~ Da, simul.infr.diff, col = "violet", lty = "solid")




lines( x =c(Da.init, 10*Da.init),  y = c(Dg.init , Dg.init + 10*Da.init))


## Observed values:


#loss in alpha vs. loss in gamma (lg transformed pearson corr)   
plot(-da, -dg ,pch = 20, col=cgam , ann=F,type ="p", axes =F, log = '')
text(x = -da,y = -dg,rownames(table.div.part) , col="grey", cex =0.6, pos = 4)
points(x = -da,y = -dg, pch = 20)
box(bty="l")
f <- cor.test(da ,dg, method = "spearman")
mtext(3, text = substitute(italic(rho == est *" "*p), list(est = round(f$est,2), p=p2star(f$p.val, marginal=F))),adj = 1, cex=0.8, font = 3, line= 0.5)
mtext(3, text = 'a)',adj = 0, cex=0.8, font = 1, line= 0.5)
mtext(1,text = expression(Delta*alpha*"-richness"), line=2 )
mtext(2,text = expression(Delta*gamma*"-richness"["c"]), line=2 , las=0)


