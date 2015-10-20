### Comparing Gamma pools at different abundance classes

### Main method : Species accumulation curves with model fitting for extrapolation

est.gamma.loss <- function(com.mat = comm, plotID = realgrasslands, spnames = impsp,
                           group=natives, threshold="th.CI", plot.graph=FALSE) {

   if (plot.graph) {
par(mfrow=c(3,4), mar=c(1,1,1,1), oma=c(2,2,1,1), las=1)
}

community <- com.mat[which((rownames(com.mat) %in% plotID)  ),]
community <- community[,colSums(community)>0]
community<- community[,colnames(community) %in% group]
community<- ceiling(community>0)


list.gammas= list()

for (sp in spnames) {
var <- com.mat[rownames(community),sp]
k <- glmSRnat.overall$impact.spread[sp, threshold]
obs <- specpool(community, var)

## use SPECACCUM to build rarefaction curve and build confidence intervals

# rare vs. at threshold:
sprare <- specaccum(community[var==1,],"random" )

# Fit a non linear model to extrapolate gamma richness below and compare to observed values above
mod <- fitspecaccum(sprare, model ="gitay" )  # fit model on observed below threshold
newdata= 1:max(obs$n)
extrapolate <-predict(mod, newdata)


### output of function
list.gammas[[sp]]<-cbind(obs,
                               data.frame(mean.null =rowMeans(extrapolate[obs$n,]) ,
                                          sd.null = apply(extrapolate[obs$n,], 1, sd),
                                          q2.5 = apply(extrapolate[obs$n,], 1, quantile, probs = 0.025),
                                          q97.5 = apply(extrapolate[obs$n,], 1, quantile, probs = 0.975),
                                          quantile = sapply(1:length(obs$n), FUN=function(x){
                                            ((sum( extrapolate[obs$n[x],]<obs$Species[x]) + sum(extrapolate[obs$n[x],] == obs$Species[x])/2) /dim(extrapolate)[2])
                                          })
                         )
)

## graphical representation
if (plot.graph) {
# above vs. below:

sp1 <- specaccum(community[var<k & var>0,],"random" )
sp2 <- specaccum(community[var>=k,],"random")
mod1 <- fitspecaccum(sp1, model ="gitay" )  # fit model on observed below threshold
extrapolate1 <-predict(mod1, newdata)
mod2 <- fitspecaccum(sp2, model ="gitay" )  # fit model on observed below threshold
extrapolate2 <-predict(mod2, newdata)

plot(sp1, ci.type="poly", col="grey", lwd=2, ci.lty=0, ci.col="lightblue", ylim=c(0,300), xlim=c(0,200))
# boxplot(sp1, col="darkgrey", add=TRUE, pch="+", border="grey")

plot.xy(xy.coords(x=newdata, y=rowMeans(extrapolate)), type = "l", lty="dashed", col="blue")
plot.xy(xy.coords(x=newdata, y=apply(extrapolate, 1, quantile,probs=0.025 )), type = "l", lty="dotted", col="blue")
plot.xy(xy.coords(x=newdata, y=apply(extrapolate, 1, quantile,probs=0.975 )), type = "l", lty="dotted", col="blue")

plot.xy(xy.coords(x=newdata, y=rowMeans(extrapolate2)), type = "l", lty="dashed", col="salmon2")
plot.xy(xy.coords(x=newdata, y=apply(extrapolate2, 1, quantile,probs=0.025 )), type = "l", lty="dotted", col="salmon2")
plot.xy(xy.coords(x=newdata, y=apply(extrapolate2, 1, quantile,probs=0.975 )), type = "l", lty="dotted", col="salmon2")
#
# plot(sp2, ci.type="poly", col="red", lwd=2, ci.lty=0, ci.col="#FF000080", ylim=c(0,300), xlim=c(0,200), add=T)
plot.xy(xy.coords(obs$n, obs$S),type="p", pch = "*")
plot.xy(xy.coords(obs$n[(k+1):length(obs$n)], obs$Species[(k+1):length(obs$n)]),
        type="p",  pch = "*",cex=2, col="red")
}
}
return(list.gammas)
}

### application for native richness
gamma.loss <- est.gamma.loss(plot.graph=FALSE)

# plot quantile results:
x11()
par(mfrow=c(3,4), mar=c(1,1,1,1), oma=c(2,2,1,1), las=1)
for (sp in spnames) {
  plot(gamma.loss[[sp]]$quantile, ylim=c(0,1))
  abline(h=0.5, lty="dotted")
}


# plot quantile results:
x11()
par(mfrow=c(3,4), mar=c(1,1,1,1), oma=c(2,2,1,1), las=1)
for (sp in spnames) {
  tmp <- gamma.loss[[sp]]
  tmp <- tmp[-1,]
  tmp <- tmp[tmp$n>4,]
  plot(as.numeric(rownames(tmp)),(tmp$Species - tmp$mean.null)/tmp$sd.null, type= "b", ylim=c(-15,8) )
  abline(h=0.5, lty="dotted")
}

x11()
par(mfrow=c(3,4), mar=c(2,2,2,1), oma=c(2,3,1,1), las=1)
for (sp in spnames) {
  tmp <- gamma.loss[[sp]][-1,]
  tmp <- tmp[tmp$n>5,]
  es <- (tmp$Species - tmp$mean.null)
  p<- c("white", "black")[as.numeric(tmp$Species < tmp$q2.5) +1]
  plot(as.numeric(rownames(tmp)),es, type= "b", xlim=c(0.5, 6.5),
       ylim=c(-max(abs(es)), max(abs(es))),pch= 21, bg = p)
  abline(h=0.5, lty="dotted")
  title(main=sp)
}
mtext(2, text = "Gamma loss vs. accum.curve(rare)", las=0, outer=T, line=1)

# ## accumulation curve of natives for all plots where species is present:
# sp1 <- specaccum(community)
# sp2 <- specaccum(community,"random" )
# plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
# boxplot(sp2, col="yellow", add=TRUE, pch="+")
# points(obs$n, obs$Species, pch = "*")
# points(obs$n[(k+1):length(obs$n)], obs$Species[(k+1):length(obs$n)], pch = "+", col="red")

# obs <- sapply(0:6, FUN = function(k){
#   if(sum(var==k)>0)  specpool(community[var<=k,]) else rep(NA, 9)
# }))
# obs.cum <- gamma.above.nat$gamma.above[sp,]

## POOLACCUM = Accumulation model => should be able to compare and even predict accumulation asymptot ?
#= = similar to specaccum, but estimate extrapolated richness indices of specpool
# in addition to number of species for random ordering of sampling units = equivalent to the random fitting nethod ?
# pool <- poolaccum(community[var<4 & var>0,], permutations = 100)
# sum.pool <- summary(pool)
# plot(pool)
#
# plot(x =sum.pool$chao[,1], y =sum.pool$S[,4], type = "l" )
# lines(x =sum.pool$S[,1], y =sum.pool$S[,2])
# lines(x =sum.pool$S[,1], y =sum.pool$S[,3])
# points(obs$n, obs$S, pch = "*")
# points(obs$n[(k+1):length(obs$n)], obs$Species[(k+1):length(obs$n)], pch = "+", col="red")
#
# plot(x =sum.pool$chao[,1], y =sum.pool$chao[,4], type = "l" )
# lines(x =sum.pool$chao[,1], y =sum.pool$chao[,2])
# lines(x =sum.pool$chao[,1], y =sum.pool$chao[,3])
# points(obs$n, obs$chao, pch = "*")
# points(obs$n[(k+1):length(obs$n)], obs$chao[(k+1):length(obs$n)], pch = "+", col="red")
#
# plot(x =sum.pool$boot[,1], y =sum.pool$boot[,4], type = "l" )
# lines(x =sum.pool$boot[,1], y =sum.pool$boot[,2])
# lines(x =sum.pool$boot[,1], y =sum.pool$boot[,3])
# points(obs$n, obs$boot, pch = "*")
# points(obs$n[(k+1):length(obs$n)], obs$boot[(k+1):length(obs$n)], pch = "+", col="red")

# estimating gamma pools sith Chao estimators

chao.gamma.nat <- function(spnames = impsp, group = natives) {
  require(vegan)

  list.gammas = list()
  for (sp in spnames) {

    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
    community <- community[,colSums(community)>0]
    community<- community[,colnames(community) %in% group]
    community<- ceiling(community>0)

    var <- comm[rownames(community),sp]

    list.gammas[[sp]] <- specpool(community, var)
  }
  return(list.gammas)
}

