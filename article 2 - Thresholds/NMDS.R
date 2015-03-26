### NMDS analysis of the invaded communities

D <- glmSR.grass$thresh

# Species which show a threshold
impsp <- rownames(D[!is.na(D$th),])

i <- 1
sp <- impsp[i]
# select realgrassland plots where at least one of the important target species is present:
community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),impsp] 
community <- comm[which((rownames(comm) %in% realgrasslands) & (rowSums(comm[,impsp])>0) ),impsp] 

nmds.results <- metaMDS(community, distance = "bray", k = 2, trymax = 10)

lim=c(min(nmds.results$points,nmds.results$species)+0.2, max(nmds.results$points,nmds.results$species)+0.2)
plot(nmds.results$points, col="darkgrey",ylim=lim, xlim=lim )
x <- nmds.results$species[,1]
y <- nmds.results$species[,2]
arrows(rep(0,length(impsp)), rep(0,length(impsp)), x, y, length =0.1, col="black")
text(x+sign(y)*0.1,y +sign(y)*0.1, label = impsp, cex=0.7 )
arrows(rep(0), rep(0), x[i], y[i], length =0.1, col="red", cex=2)
