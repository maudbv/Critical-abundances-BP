### trying to enerate null distributions of species

# characterizing the distribution of observations per abudnance classes

tmp=data.frame(as.matrix(table(databp$SpeciesCode, databp$abun))[,1:6]) # nummber of obs in each ab classes per species
tmp$zero=as.vector(apply(comm, 2, function(x) sum(x==0))[rownames(tmp)]) # number of absences for each species
tmp=tmp[,c(7,1:6)] # reorder matrix

m=apply(tmp,2,mean) # mean frequencies per class
e=apply(tmp,2,function(x) sd(x)/sqrt(length(x))) # standard error of frequencies
plot(1:7,m,log="y", type="h", lwd=15, col="grey", lend=1)
segments(1:7, m-e, 1:7, m+e )


hist(tmp[,1])
hist(tmp[,2])
hist(tmp[,3])
hist(tmp[,4])
hist(tmp[,5])
hist(tmp[,6])
hist(tmp[,7])

dnorm(0:200,m[2], e[2])
