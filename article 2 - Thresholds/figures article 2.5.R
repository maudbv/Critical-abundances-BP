# Define some elements
abclasses= c("Rare" ,"Occasional",  "Frequent", "Common", "Abundant","Dominant")
threshold ="th.CI"
sel <- impsp


### Figure 1: barplot of frequencies:
effects =list(glmSRnat.sum$class.summary,glmSRali.sum$class.summary )

n <- length(unique(effects[[1]]$group))

par(mfcol=c(n,2), oma=c(5,6,4,1), mar=c(1,1,1,2), las=1)
for (j in 1:2) {
  sum.df = effects[[j]]
  for (i in 1:n) {
    S <- as.data.frame(sum.df[sum.df$group == unique(sum.df$group)[i],])
    barplot(S$nb.sp, ylim=c(0,max(40, S$nb.sp)),col= "grey",  border= NA, axes=F)

    par(new=T)
    barplot(S$freq.negative.above, ylim=c(0,max(40, S$nb.sp)),col= "grey50",  border= NA, axes=F)

    par(new=T)
    b <- barplot(S$freq.thr, col="black" , ylim=c(0,max(40, S$nb.sp)), border= NA, axes=F)

    axis(2, tcl = -0.3, cex.axis = 0.8, mgp = c(1,0.5,0), las = 1)
    if (i==2) text(y=-3, x = b, labels= abclasses[2:6],  cex=0.8, srt=45, adj=1, xpd = NA)

    if (j==2 & i==1) {
      legend(x=0.3, y=40, bty="n", bg="white",legend=c("Critical abundance", "Total occurrences"),
             fill=c("black",  "grey"), border= c("black",  "grey"), cex=0.9, xpd=NA)
    }
    if(j==1) {
      mtext(text="number of species", side=2, outer=F, line=1.5, las=0, cex=0.85)
    }
  }
}

mtext(side=2, text=c("Alien\nfocal\nspecies", "Native\nfocal\nspecies"),line=4, at=c(0.25,0.75),adj=0.5, outer=T, las=1)
mtext(side=3, text=c("a) Native richness", "b) Alien richness"),line=0, at=c(0.25,0.75),adj=0.5, outer=T, las=1)

mtext(text="Abundance class", side=1, outer=T, line=3)