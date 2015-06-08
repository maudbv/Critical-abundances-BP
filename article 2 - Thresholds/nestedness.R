### nestedness and betadiversity of communities
library(vegan)
 
## Overall partitionning of diversity
div.part.overall <- ( function ()
  {
community <- comm[which(rownames(comm) %in% realgrasslands),]
community <- community[,colSums(community)>0]
community<- ceiling(community>0)

n <- sum(rowSums(community)>0)
gam <- sum(colSums(community)>0)
alpha <- mean(rowSums(community), na.rm=T)
beta.add <- gam - alpha
beta.mult <- gam/alpha
beta.prop <- (gam - alpha)/gam

div.part.overall <- c(n, gam, alpha, beta.add, beta.mult,beta.prop)

comm.nat <- community[,colnames(community) %in% natives] 

n <- sum(rowSums(comm.nat)>0)
gam <- sum(colSums(comm.nat)>0)
alpha <- mean(rowSums(comm.nat), na.rm=T)
beta.add <- gam - alpha
beta.mult <- gam/alpha
beta.prop <- (gam - alpha)/gam

div.part.overall <- rbind(div.part.overall,
                          c(n, gam, alpha, beta.add, beta.mult,beta.prop))

comm.ali <- community[,colnames(community) %in% aliens] 
n <- sum(rowSums(comm.ali)>0)
gam <- sum(colSums(comm.ali)>0)
alpha <- mean(rowSums(comm.ali), na.rm=T)
beta.add <- gam - alpha
beta.mult <- gam/alpha
beta.prop <- (gam - alpha)/gam

div.part.overall <- rbind(div.part.overall,
                          c(n, gam, alpha, beta.add, beta.mult,beta.prop))
colnames(div.part.overall) = c("nb.plot", "gamm", "alpha", "beta.add", "beta.mult","beta.prop")
rownames(div.part.overall) = c("total richness", "native richness", "alien richness")
 return(div.part.overall)
})()


# Partition Gamma/beta/alpha native diversity for each imp species above/below

div.part.nat <- matrix(NA, nrow = 11, ncol = 6*2,
                   dimnames= list(
                     impsp,
                     paste(c(rep("below",6), rep("above", 6)),
                           c("nb.plot", "gamm", "alpha", "beta.add", "beta.mult","beta.prop"),
                           sep=".")))

for (sp in impsp) {
  
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,colSums(community)>0]
  community<- community[,colnames(community) %in% natives] 
  community<- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  n <- sum(rowSums(community[var==0,])>0)
  gam <- sum(colSums(community[var==0,])>0)
  alpha <- mean(rowSums(community[var==0,]), na.rm=T)
  beta.add <- gam - alpha
  beta.mult <- gam/alpha
  beta.prop <- (gam - alpha)/gam
    
  n2 <- sum(rowSums(community[var==1,])>0)
  gam2 <- sum(colSums(community[var==1,])>0)
  alpha2 <- mean(rowSums(community[var==1,]), na.rm=T)
  beta.add2 <- gam2 - alpha2
  beta.mult2 <- gam2/alpha2
  beta.prop2 <- (gam2 - alpha2)/gam2
  
  div.part.nat[sp,] <- c(n, gam, alpha, beta.add, beta.mult,beta.prop,
                    n2, gam2, alpha2, beta.add2, beta.mult2,beta.prop2)
  
}

# Partition Gamma/beta/alpha ALIEN diversity for each imp species above/below

div.part.ali <- matrix(NA, nrow = 11, ncol = 6*2,
                       dimnames= list(
                         impsp,
                         paste(c(rep("below",6), rep("above", 6)),
                               c("nb.plot", "gamm", "alpha", "beta.add", "beta.mult","beta.prop"),
                               sep=".")))

for (sp in impsp) {
  
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,colSums(community)>0]
  community<- community[,colnames(community) %in% aliens] 
  community<- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  n <- sum(rowSums(community[var==0,])>0)
  gam <- sum(colSums(community[var==0,])>0)
  alpha <- mean(rowSums(community[var==0,]), na.rm=T)
  beta.add <- gam - alpha
  beta.mult <- gam/alpha
  beta.prop <- (gam - alpha)/gam
  
  n2 <- sum(rowSums(community[var==1,])>0)
  gam2 <- sum(colSums(community[var==1,])>0)
  alpha2 <- mean(rowSums(community[var==1,]), na.rm=T)
  beta.add2 <- gam2 - alpha2
  beta.mult2 <- gam2/alpha2
  beta.prop2 <- (gam2 - alpha2)/gam2
  
  div.part.ali[sp,] <- c(n, gam, alpha, beta.add, beta.mult,beta.prop,
                         n2, gam2, alpha2, beta.add2, beta.mult2,beta.prop2)
  
}

#compare species list  

## check that mean frequency of species is comparable between above and below
i=2
sp=impsp[i]
community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
community <- ceiling(community>0)
community <- community[,colSums(community)>0]


comm.nat<- community[,colnames(community) %in% natives] 
comm.ali<- community[,colnames(community) %in% aliens] 

var <- comm[rownames(community),sp]
var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1

cdat <- as.matrix(comm.nat[var==0,])
indices <-seq(10,length(rownames(cdat)),10)

m <- s <- NULL
for (k in indices){
tmp <- sapply(1:10, function(x){
       samp <- colSums(cdat [sample(x = rownames(cdat) ,k, replace =F),])
       return(samp)
})
m <- rbind(m,rowMeans(tmp))
s <- rbind(s, apply(tmp, 1, sd))
}
  plot(indices, seq(0,180, length.out = length(indices)), type= "n" )     
for(i in colnames(m)){
  points(indices, m[,i], type= "b")    
  segments(indices, m[,i] -s[,i], indices, m[,i] +s[,i], col="grey")
}
cor.test(m[1,], m[57,], method = "spearman")
plot(m[1,]/10, m[57,]/570)
abline(0,1)

# Compare frequency of natives above and below

#community change plots
## species order
cm <- comm[(rownames(comm) %in% realgrasslands),]
cm <- cm[,colnames(cm)%in%natives]
cm <- ceiling(cm>0)
cm <- cm[,colSums(cm)>0]
ord <- colnames(cm)[order(colSums(cm))]
ord <- ord[ord %in% natives]

x11()
par( mar=c(1,3,1,0))
y =  colSums(cm)[order(-colSums(cm))]/length(rowSums(cm)>0)
names(y) <- species[names(y),]$tip

plot(1:length(y), y,axes= F, ann=F, type= "h", col="grey")
text((1:12)+2, y[1:12],labels= names(y)[1:12], adj=0, cex=0.7)
axis(2, tcl = -0.3, las=1, cex = 0.7)
t(data.frame( t(y)))

x11()
 par(mfrow=c(3,4), mar=c(2,1,3,1))

for (sp in impsp) {

  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,colSums(community)>0]
  community<- community[,colnames(community) %in% natives] 
  
  comm.nat <- community
  comm.nat <- ceiling(comm.nat>0)
  
var <- comm[rownames(community),sp]
var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1

B.freq <- colSums(comm.nat[var==0,])/sum(var==0)
A.freq <- colSums(comm.nat[var==1,])/sum(var==1)


cvec <-  c("grey","red")[(A.freq == 0) + 1]
sum((A.freq == 0) )
cvec [B.freq == 0] <- "forestgreen"
sum(B.freq == 0) 
names(cvec) <- names(B.freq)

f1 <- wilcox.test(  A.freq ,B.freq, paired=T, alternative = "less")

## community change plots

#unordered logs
# ord <- order(colSums(comm.nat))
    
barplot(log(0.005)-log(B.freq[ord]+0.005), las=2, xlim=c(-7,7),
        cex.names = 0.5, horiz = T, col = cvec[ord],
        border = cvec[ord], xaxt= "n", yaxt = "n")

barplot(-(log(0.005)-log(A.freq[ord]+0.005)), las=2, xlim=c(-7,7),
        cex.names = 0.5, horiz = T, col = cvec[ord], 
        border = cvec[ord], add=T, xaxt= "n", yaxt ="n")

mtext(2, text = sum((A.freq == 0) ), line = -1, las=1, col= "red", cex= 0.7)

mtext(4, text = sum((B.freq == 0) ), line = -1, las=1, col= "forestgreen", cex= 0.7)

mtext(1, text = paste(sum((B.freq!= 0) & A.freq!= 0), "shared"), line = -2, adj=0.1,las=1, col= "black", cex= 0.7)

# #logs
# barplot(log(0.005)-log(sort(B.freq+0.005)), las=2, xlim=c(-7,7),
#         cex.names = 0.5, horiz = T, col = cvec[order(B.freq)],
#         border = cvec[order(B.freq)], xaxt= "n", yaxt = "n")
# 
# barplot(-(log(0.005)-log(A.freq[order(B.freq)]+0.005)), las=2, xlim=c(-7,7),
#         cex.names = 0.5, horiz = T, col = cvec[order(B.freq)], 
#         border = cvec[order(B.freq)], add=T, xaxt= "n", yaxt ="n")
# 
# #no logs
# barplot(-sort(B.freq), las=2, xlim=c(-0.5,0.5),
#         cex.names = 0.5, horiz = T, col = cvec[order(B.freq)],
#         border = cvec[order(B.freq)], xaxt= "n", yaxt = "n")
# 
# barplot(A.freq[order(B.freq)], las=2, xlim=c(-0.5,+0.5),
#         cex.names = 0.5, horiz = T, col = cvec[order(B.freq)], 
#         border = cvec[order(B.freq)], add=T, xaxt= "n", yaxt ="n")
abline(v=0)

#axis(1, at = c(-0.4,-0.20,0,0.2,0.4),label= c("40%","20%",0,"20%","40%"), tcl = 0.2,mgp=c(1,0.2,0))
axis(1, at = (log(0.005)-log(c(0.7,0.1,0.01,0.01,0.1,0.7)+0.005))*c(1,1,1,-1,-1,-1),
     label= c("70%","10%","1%","1%","10%","70%"), tcl = 0.2, 
     mgp=c(1,0.2,0))

mtext(side = 3, text =c("below", "above"),adj = c(0.33,+0.66), cex = 0.7, font = 3)

mtext(sp, line = 1.2, cex = 0.7)

}


### ALIEN POOLS
cm <- comm[(rownames(comm) %in% realgrasslands),]
cm <- cm[,colnames(cm)%in%aliens]
cm <- ceiling(cm>0)
cm <- cm[,colSums(cm)>0]
ord <- colnames(cm)[order(colSums(cm))]
ord <- ord[ord %in% aliens]

x11()
par( mar=c(1,3,1,0))
y =  colSums(cm)[order(-colSums(cm))]/length(rowSums(cm)>0)
names(y) <- species[names(y),]$tip

plot(1:length(y), y,axes= F, ann=F, type= "h", col="grey")
text((1:22)+2, y[1:22],labels= names(y)[1:22], adj=0, cex=0.7)
axis(2, tcl = -0.3, las=1, cex = 0.7)
t(data.frame( t(y)))

x11()
par(mfrow=c(3,4), mar=c(2,1,3,1))

for (sp in impsp) {
  
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,colSums(community)>0]
  community<- community[,colnames(community) %in% aliens] 
  
  comm.ali <- community
  comm.ali <- community[,-grep(sp, names(comm))]
  comm.ali <- ceiling(comm.ali>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  B.freq <- colSums(comm.ali[var==0,])/sum(var==0)
  A.freq <- colSums(comm.ali[var==1,])/sum(var==1)
  
  
  cvec <-  c("grey","red")[(A.freq == 0) + 1]
  sum((A.freq == 0) )
  cvec [B.freq == 0] <- "forestgreen"
  sum(B.freq == 0) 
  names(cvec) <- names(B.freq)
  
  f1 <- wilcox.test(  A.freq ,B.freq, paired=T, alternative = "less")
  
  ## community change plots
  
  #unordered logs
  # ord <- order(colSums(comm.nat))
  
  barplot(log(0.005)-log(B.freq[ord]+0.005), las=2, xlim=c(-7,7),
          cex.names = 0.5, horiz = T, col = cvec[ord],
          border = cvec[ord], xaxt= "n", yaxt = "n")
  
  barplot(-(log(0.005)-log(A.freq[ord]+0.005)), las=2, xlim=c(-7,7),
          cex.names = 0.5, horiz = T, col = cvec[ord], 
          border = cvec[ord], add=T, xaxt= "n", yaxt ="n")
  
  abline(v=0)
  
  mtext(2, text = sum((A.freq == 0) ), line = -1, las=1, col= "red", cex= 0.7)
  
  mtext(4, text = sum((B.freq == 0) ), line = -1, las=1, col= "forestgreen", cex= 0.7)
  
  mtext(1, text = paste(sum((B.freq!= 0) & A.freq!= 0), "shared"), line = -2, adj=0.1,las=1, col= "black", cex= 0.7)
  
  #axis(1, at = c(-0.4,-0.20,0,0.2,0.4),label= c("40%","20%",0,"20%","40%"), tcl = 0.2,mgp=c(1,0.2,0))
  axis(1, at = (log(0.005)-log(c(0.7,0.1,0.01,0.01,0.1,0.7)+0.005))*c(1,1,1,-1,-1,-1),
       label= c("70%","10%","1%","1%","10%","70%"), tcl = 0.2, 
       mgp=c(1,0.2,0))
  
  mtext(side = 3, text =c("below", "above"),adj = c(0.33,+0.66), cex = 0.7, font = 3)
  
  mtext(sp, line = 1.2, cex = 0.7)
  
}

# correlation graphs
x11()
par(mfrow=c(3,4), mar=c(2,2,2,2), oma= c(2,3,1,1))
for (sp in impsp) {
  
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- ceiling(community>0)
  community <- community[,colSums(community)>0]
  
  comm.nat<- community[,colnames(community) %in% natives] 
  comm.ali<- community[,colnames(community) %in% aliens] 
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  B.freq <- colSums(comm.nat[var==0,])/sum(var==0)
  A.freq <- colSums(comm.nat[var==1,])/sum(var==1)
   
  cvec <-  c("grey","red")[(A.freq == 0) + 1]
  sum((A.freq == 0) )
  cvec [B.freq == 0] <- "forestgreen"
  sum(B.freq == 0) 
  
  f1 <- wilcox.test(  A.freq ,B.freq, paired=T, alternative = "less")
  c1 <- cor.test( A.freq , B.freq, method = "spearman")
  
  #correlation plot
#   plot(jitter(B.freq, factor = 1), jitter(A.freq, factor=1),
#        ylim=c(-0.01,0.9), xlim=c(-0.01,0.9),pch=21, bg =cvec,
#        ann=F, axes =F)
  plot(B.freq+0.005, A.freq+0.005,log="xy",pch=21, bg =cvec, ann=F, axes =F,
       ylim=c(0.005, 0.9), xlim=c(0.005,0.9))
  
  abline(0,1)
  mtext(paste(sp, "\n rho=", round(c1$est, 2), p2star(c1$p.val)), cex = 0.7)

#   axis(1, at =c(0.005,0.01,seq(0.1,0.9,0.1)),label=c(0,0.01,seq(10,90,10)), tcl = 0.2,mgp=c(1,0.2,0), cex=0.7)
#   axis(2, at =c(0.005,0.01,seq(0.1,0.9,0.1)),label=c(0,0.01,seq(10,90,10)), tcl = 0.2,mgp=c(1,0.2,0),las=1,  cex=0.7)
  axis(1, tcl = 0.2,mgp=c(1,0.2,0), cex=0.7)
 axis(2,tcl = 0.2,mgp=c(1,0.2,0),las=1,  cex=0.7)

  box(bty="l")    
}
mtext(side = 1, text =c("below critical"),cex = 1, line=0,outer= T)
mtext(side = 2, text =c("above critical"),cex = 1, line= 1,outer= T)


##other types of graphs
x11()
par(mfrow=c(3,4), mar=c(2,1,3,1))
for (sp in impsp) {
  
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- ceiling(community>0)
  community <- community[,colSums(community)>0]
  
  comm.nat<- community[,colnames(community) %in% natives] 
  comm.ali<- community[,colnames(community) %in% aliens] 
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  B.freq <- colSums(comm.nat[var==0,])/sum(var==0)
  A.freq <- colSums(comm.nat[var==1,])/sum(var==1)
  
  cvec <-  c("grey","red")[(A.freq == 0) + 1]
  sum((A.freq == 0) )
  cvec [B.freq == 0] <- "forestgreen"
  sum(B.freq == 0) 
  
  f1 <- wilcox.test(  A.freq ,B.freq, paired=T, alternative = "less")
  c1 <- cor.test( A.freq , B.freq, method = "spearman")
  
  # correlation plot
  # plot(jitter(B.freq, factor = 2), A.freq, ylim=c(0,0.5), xlim=c(-0.01,0.5),pch=21, bg =cvec)
  # abline(0,1)
  # mtext(paste(sp, "\n rho=", round(c1$est, 2), p2star(c1$p.val)))
  
  # # deviation plot
  # deltafreq <- (A.freq - B.freq)
  # barplot(sort(deltafreq), las=2,
  #         cex.names = 0.5, horiz = T, col = cvec[order(deltafreq)], 
  #         border = cvec[order(deltafreq)])
  # mtext(paste(sp, "\n decrease",p2star(f1$p.val)))
  
  
  ## community plots
  
  barplot(-sort(B.freq), las=2, xlim=c(-0.5,0.5),
          cex.names = 0.5, horiz = T, col = cvec[order(B.freq)],
          border = cvec[order(B.freq)], xaxt= "n", yaxt = "n")
  
  barplot(A.freq[order(B.freq)], las=2, xlim=c(-0.5,+0.5),
          cex.names = 0.5, horiz = T, col = cvec[order(B.freq)], 
          border = cvec[order(B.freq)], add=T, xaxt= "n", yaxt ="n")
  abline(v=0)
  
  axis(1, at = c(-0.4,-0.20,0,0.2,0.4),label= c("40%","20%",0,"20%","40%"), tcl = 0.2, 
       mgp=c(1,0.2,0))
  mtext(side = 3, text =c("below", "above"),at = c(-0.2,+0.2), cex = 0.7, font = 3)
  
  mtext(sp, line = 1.2, cex = 0.7)
  
  # paired t tplot
  # plot(c(rep(1,length(B.freq)), rep(2, length(A.freq))),
  #        c(B.freq, A.freq),
  #        ylim=c(0,0.5),xlim=c(0.5,2.5),
  #      type= "p", pch=21, bg =cvec)
  # segments(rep(1,length(B.freq)), 
  #  B.freq, 
  #   rep(2, length(A.freq)),
  #  A.freq,
  # bg =cvec)
  
}

# natives not significantly nested
checker.SRnat <- oecosimu(comm.nat, nestedchecker, method = "r1", nsimul = 599, burnin = 0, thin = 1,
         statistic = "statistic", alternative = c("two.sided", "less", "greater"))
  
nestbetasor.SRnat <- oecosimu(comm.nat, nestedbetasor, method = "r1", nsimul = 99, burnin = 0, thin = 1,
                              statistic = "statistic", alternative = c("two.sided", "less", "greater"))

#alien coms    
nestbetasor.SRali <- oecosimu(comm.ali, nestedbetasor, method = "r1", nsimul = 99, burnin = 0, thin = 1,
                              statistic = "statistic", alternative = c("two.sided", "less", "greater"))

checker.SRali <- oecosimu(comm.ali, nestedchecker, method = "r1", nsimul = 99, burnin = 0, thin = 1,
                         statistic = "statistic", alternative = c("two.sided", "less", "greater"))


### temperature
out <- nestedtemp(comm.ali)
out
plot(out)
plot(out, kind="incid")
