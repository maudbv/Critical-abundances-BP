## TRait variation in banks peninsula

## Alien species traits

alientraits=traitdata[which(traitdata$ALIEN==1),]
alientraits$Nat_end - alientraits$Nat_start

recent=alientraits[which(alientraits$Nat_start>1900),]
alientraits=orderBy( ~ Nat_start + Nat_end , data=alientraits)

pdf("period of naturalisation.pdf")
plot(x=1, y=1, type="n", xlab=NA,las=1, xlim=c(1840,2000), ylim=c(0,280),
    main="Period of naturalisation", yaxt="n", ylab=NA)
segments(alientraits$Nat_start, 1:311,alientraits$Nat_end, 1:311, lwd=1,
         col=c("grey", "grey30")[as.numeric(rownames(alientraits) %in% rownames(effect.summary))+1])
segments(alientraits$Nat_start, 1:311,alientraits$Nat_end, 1:311, lwd=1,
         col=c(NA, "firebrick")[as.numeric(rownames(alientraits) %in% rownames(effect.summary[effect.summary$total.effect>0,]))+1])
text(x=sort(unique((alientraits$Nat_end + alientraits$Nat_start)/2)),
     y=cumsum(table(alientraits$Nat_start))+10 ,
     label=as.numeric(table(alientraits$Nat_start)))
legend("bottomright", lty=c("solid", "solid", "solid"), col=c("grey", "grey30", "firebrick"),
       legend=c("alien", "target alien", "significant alien"), bty="n")
dev.off()

is.factor(alientrait$Naturalisation)
barplot(table(alientraits[, "Naturalisation"]))
barplot(table(alientraits[rownames(alientraits) %in% rownames(effect.summary), "Naturalisation"]))
barplot(table(alientraits[, "Naturalisation"]))

# effect of naturalisation time on correlations


table(effect.summary$effect[effect.summary$ALIEN==1], effect.summary$Nat[effect.summary$ALIEN==1])

plot(S.pval ~ as.factor(Nat), data= effect.summary[effect.summary$ALIEN==1,])

# add naturalisation in effect summary table
effect.summary$Nat=alientraits[rownames(effect.summary),c("Naturalisation") ]
tmp=effect.summary[, c("SpeciesCode","Division","Class","Family","SpeciesName","Growth.forms","ALIEN","WEEDOC",
                       "Sp.occurence.x","nb.class","max","min",
                       "effect","rho","S.pval",
                       "nat.effect","nat.rho","nat.S.pval",
                       "ali.effect","ali.rho","ali.S.pval",
                       "total.effect")]

write.csv(tmp, file="effect.summary.csv")

## add in species
species$Nat=alientraits[rownames(species),c("Naturalisation") ]
species$Nat.start=alientraits[rownames(species),c("Nat_start") ]
species$Nat.end=alientraits[rownames(species),c("Nat_end") ]



