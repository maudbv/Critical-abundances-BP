# trying to find a good impact index combining impact size and prevalence

x = seq(0,1, 0.01)
S <- matrix(x, length(x),length(x), byrow=T)
P <- matrix(x, length(x),length(x),byrow=F)

z1 <- (S+P)/sqrt(2)
z2 <-  S*P
z3 <-  sqrt(S*P)

# 
# plot(1:(length(x)^2),z0, col=rainbow(11))
# plot(1:(length(x)^2),z1, col=rainbow(11))
# plot(1:(length(x)^2),z2, col=rainbow(11))
# 
# plot(z1,z2, col=rainbow(11))
# plot(z2,z1, col=rainbow(11))

library(lattice)

# best might be simple multiplication :
levelplot(z2~ S * P, col.regions= colorRampPalette(c("beige" , "firebrick")), main = "I = S * P" )
levelplot(z3~ S * P, col.regions= colorRampPalette(c("beige",  "firebrick")) , main = "I = sqrt(S*P)")

