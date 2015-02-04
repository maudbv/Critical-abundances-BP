
test.sim1rhos=function(nm=sim1.rhos, subsample="none", graph=F) {

if(subsample=="A") nm=nm[which(rownames(nm) %in% rownames(species)[which(species$ALIEN==1)]),]
if(subsample=="N") nm=nm[which(rownames(nm) %in% rownames(species)[which(species$ALIEN==0)]),]

if (graph) {
  plot(density(na.omit(unlist(nm))), col="blue", main=NA )
 abline(v=mean(na.omit(unlist(nm))), col="blue")
}
nm1rhos.sig=rep(NA,dim(nm)[2])

for (i in 2:dim(nm)[2]){
x=data.frame( rho= nm[,1], null=nm[,i])
rownames(x) = rownames(nm)

if (graph) lines(density(na.omit(x$null)), col="grey")
  
y= quantile(x$rho, 0.75, na.rm=T) +0.1
t=wilcox.test(x$rho,x[,2], paired=T)
nm1rhos.sig[i]=t$p.val

print(i)
# # 
# if (ylab) mtext(side=2, text=expression("Spearman's "*rho), outer=F, las=0, line=1.2, cex=0.65)
}

if (graph) {
  lines(density(na.omit(nm[,1])), col="red")
abline(v=mean(na.omit(nm[,1])), col="red")
}

out=sum(nm1rhos.sig>0.05, na.rm=T)/dim(nm)[2]
print(out)
return(out)
}


###### run the function ############
rhotest=matrix(rep(NA, 27), nrow=3, ncol=9)
rownames(rhotest)=c("all", "grass","wood")
colnames(rhotest)=c("SR", "NR", "AR", "SRnat", "NRnat", "ARnat", "SRali", "NRali", "ARali")

x11()
par(mfrow=c(3,3), oma=c(0,0,2,0))

for (i in 1:3) {
  subsample=c("none", "N", "A")[i]
  for (j in 1:3) {
    nm=list(sim1.rhos,sim1.rhos.nat, sim1.rhos.ali)[[j]]
    rhotest[1, (i-1)*3 + j]=test.sim1rhos(nm, subsample, graph=T)
    title(main=paste(subsample, p2star(rhotest[1, (i-1)*3 + j]))) 
  }
}
mtext(side=3, text="All plots", outer=T)


x11()
par(mfrow=c(3,3), oma=c(0,0,2,0))

for (i in 1:3) {
  subsample=c("none", "N", "A")[i]
  for (j in 1:3) {
    nm=list(sim1.rhos.grass,sim1.rhos.nat.grass, sim1.rhos.ali.grass)[[j]]
    rhotest[2, (i-1)*3 + j]=test.sim1rhos(nm, subsample, graph=T)
    title(main=paste(subsample, p2star(rhotest[2, (i-1)*3 + j]))) 
  }
}
mtext(side=3, text="Grasslands", outer=T)

x11()
par(mfrow=c(3,3), oma=c(0,0,2,0))
for (i in 1:3) {
  subsample=c("none", "N", "A")[i]
  for (j in 1:3) {
    nm=list(sim1.rhos.wood,sim1.rhos.nat.wood, sim1.rhos.ali.wood)[[j]]
    rhotest[3, (i-1)*3 + j]=test.sim1rhos(nm, subsample, graph=T)
    title(main=paste(subsample, p2star(rhotest[3, (i-1)*3 + j]))) 
  }
}
mtext(side=3, text="Woodlands", outer=T)

### illustrate why other nthod was wrong
nm=sim1.rhos.ali.wood
nm=nm[rownames(nm) %in% aliens,]
test.sim1rhos(nm=nm, subsample="A", graph=T)
abline(v=mean(unlist(nm[rownames(nm) %in% aliens,]), na.rm=T))
abline(v=mean(nm[rownames(nm) %in% aliens,1], na.rm=T), col="red")
lines(density(rowMeans(nm, na.rm=T)), col="blue")


plot(density(na.omit(nm[,1])), col="red", ylim=c(0,15))
lines(density(na.omit(unlist(nm))))
lines(density(rowMeans(nm, na.rm=T)), col="blue")


plot(rep(2,dim(nm)[1]),nm[,1], xlim=c(0.5,2.5), type="p")
abline(h=0, col="grey")
points(rep(1,dim(nm)[1]) , rowMeans(nm, na.rm=T))
segments(rep(1,dim(nm)[1]) , rowMeans(nm, na.rm=T), rep(2,dim(nm)[1]),nm[,1], col="darkgrey")


### Null distribution for all species
x11()
par(mfrow=c(3,3), oma=c(0,0,2,0))

for (i in 1:3) {
  subsample=c("none", "N", "A")[i]
  for (j in 1:3) {
    nm=list(sim1.rhos,sim1.rhos.nat, sim1.rhos.ali)[[j]]
    
    plot(density(as.numeric(nm[1,])), na.rm=T, ylim=c(0,6))
    for (i in 2:dim(nm)[1]) {
      lines(density(as.numeric(nm[i,]), na.rm=T))
    }
    
    title(main=paste(subsample) )
  }
}
mtext(side=3, text="All plots", outer=T)


######### t tests
rhot.test= matrix(rep(NA, 27), nrow=3, ncol=9)
rownames(rhot.test)=c("all", "grass","wood")
colnames(rhot.test)=c("SR", "NR", "AR", "SRnat", "NRnat", "ARnat", "SRali", "NRali", "ARali")

for (i in 1:3) {
  subsample=c("none", "N", "A")[i]
  for (j in 1:3) {
    nm=list(sim1.rhos,sim1.rhos.nat, sim1.rhos.ali)[[j]]
    if(subsample=="A") nm=nm[which(rownames(nm) %in% rownames(species)[which(species$ALIEN==1)]),]
    if(subsample=="N") nm=nm[which(rownames(nm) %in% rownames(species)[which(species$ALIEN==0)]),]
    rhot.test[1, (i-1)*3 + j]=round(t.test(nm[,1])$p.val, 4)
    rhot.test[1, (i-1)*3 + j]=round(wilcox.test(nm[,1])$p.val, 4)
  }
}

for (i in 1:3) {
  subsample=c("none", "N", "A")[i]
  for (j in 1:3) {
    nm=list(sim1.rhos.grass,sim1.rhos.nat.grass, sim1.rhos.ali.grass)[[j]]
    if(subsample=="A") nm=nm[which(rownames(nm) %in% rownames(species)[which(species$ALIEN==1)]),]
    if(subsample=="N") nm=nm[which(rownames(nm) %in% rownames(species)[which(species$ALIEN==0)]),]
    rhot.test[2, (i-1)*3 + j]=round(t.test(nm[,1])$p.val, 4)
    rhot.test[2, (i-1)*3 + j]=round(wilcox.test(nm[,1])$p.val, 4)
  }
}

for (i in 1:3) {
  subsample=c("none", "N", "A")[i]
  for (j in 1:3) {
    nm=list(sim1.rhos.wood,sim1.rhos.nat.wood, sim1.rhos.ali.wood)[[j]]
    if(subsample=="A") nm=nm[which(rownames(nm) %in% rownames(species)[which(species$ALIEN==1)]),]
    if(subsample=="N") nm=nm[which(rownames(nm) %in% rownames(species)[which(species$ALIEN==0)]),]
    rhot.test[3, (i-1)*3 + j]=round(t.test(nm[,1])$p.val, 4)
    rhot.test[3, (i-1)*3 + j]=round(wilcox.test(nm[,1])$p.val, 4)
  }
}

