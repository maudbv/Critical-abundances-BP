### representing NM1 and NM2 tests on frequency of significant relationships

nm=sim1

plot(table(nm$significant), type="h", ylim=c(0,400))
points(1:length(table(nm$positive))+ 0.1, table(nm$positive), type="h", col="forestgreen")
points(1:length(table(nm$negative)) - 0.2, table(nm$negative), type="h", col="red")
points(nm$significant[1],50,pch="*")
points(nm$positive[1], pch="*" col="forestgreen")
points(1:length(table(nm$negative)) - 0.2, table(nm$negative), type="h", col="red")

nm=shuff.all

plot(table(nm$Natives.negSR),  ylim=c(-10,400), type="h",lwd=5, lend=1, 
     xlab="number of species with significant positive relationships",
     ylab="Null frequency" )
points(1:length(table(nm$Aliens.negSR))+ 0.1, table(nm$Aliens.negSR),
       type="h", lwd=5, lend=1, col="forestgreen")
points(nm$Natives.negSR[1],-10,pch=24)
points(nm$Aliens.negSR[1],-10,pch=24, col="forestgreen", bg="forestgreen")

plot(table(nm$Natives.posSR),  type="h", lwd=5, lend=1, 
     xlab="number of species with significant negative relationships",
     ylab="Null frequency" ,
     ylim=c(0,400), xlim=c(0,10))
points(1:length(table(nm$Aliens.posSR))+ 0.1, table(nm$Aliens.posSR), 
       type="h", lwd=5, lend=1, col="forestgreen")

points(nm$Aliens.posSR[1],-10,pch=24,col="forestgreen", bg="forestgreen")
points(nm$Natives.posSR[1],-10,pch=24)


## NAtive richness
x11()

plot(table(nm$Natives.negSRnat),  ylim=c(-10,300), type="h",lwd=3, lend=1, 
     main="Significant negative relationships", xlab="number of species",
     ylab="Null frequency" )
points(1:length(table(nm$Aliens.negSRnat))+ 0.1, table(nm$Aliens.negSRnat),
       type="h", lwd=5, lend=1, col="forestgreen")
points(nm$Natives.negSRnat[1],-10,pch=24)
points(nm$Aliens.negSRnat[1],-10,pch=24, col="forestgreen", bg="forestgreen")

legend(8,300,, legend = c("Native targets", "Alien targets"), 
        fill=c("black", "forestgreen"), bty="n", cex=0.8)
legend(8.1,255, legend = "observed value", pch=24, bty="n", cex=0.8)

plot(table(nm$Natives.posSRnat),  type="h", lwd=4, lend=1, 
     main="Significant positive relationships",xlab="number of species",
     ylab="Null frequency",
     ylim=c(-10,400))
points(1:length(table(nm$Aliens.posSRnat))+ 0.1, table(nm$Aliens.posSRnat), 
       type="h", lwd=5, lend=1, col="forestgreen")
points(nm$Aliens.posSRnat[1],-10,pch=24,col="forestgreen", bg="forestgreen")
points(nm$Natives.posSRnat[1],-10,pch=24)