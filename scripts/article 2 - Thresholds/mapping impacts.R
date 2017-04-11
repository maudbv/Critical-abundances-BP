
### Mapping impacts

## proportion of invasion AR/TR
envplot$prop.alien <- envplot$SRali/envplot$SR
x <- cut(envplot$prop.alien , 20)
cp <- colorRampPalette(c("white" ,"tan","sienna4"))(20)
plot(study_area, border = "grey", col="grey")
plot(envplot, add = T, pch=22, col = cp[as.numeric(x)], bg = cp[as.numeric(x)])

## dominant species
envplot$first.rank <-NA
firstranksp <- databp[databp$DominanceRank == 1 , c("SpeciesCode", "PlotName")]
firstranksp$ALIEN <- 0
firstranksp$ALIEN[firstranksp$SpeciesCode %in% aliens] <- 1
envplot$first.rank <-firstranksp$SpeciesCode[match(envplot$PLOTID,firstranksp$PlotName)]
envplot$ALIEN.dom <-firstranksp$ALIEN[match(envplot$PLOTID,firstranksp$PlotName)]

tmp.envplot <-envplot[envplot$PLOTID %in% realgrasslands,]

plot(study_area, border = "grey", col="grey")
points(tmp.envplot,
       pch=22, col=c("goldenrod", "firebrick")[tmp.envplot$ALIEN.dom +1],
       bg =c("goldenrod", "firebrick")[tmp.envplot$ALIEN.dom +1])

### Impact spread :
pdf(file = "map impact.pdf")
db=databp[databp$PlotName %in% realgrasslands,]
targets <- impsp

par(mfrow=c(3,4), mar=c(0,0,0,0),oma= c(2,0,0,0))

for( i in targets) {
  plot.alien=as.character(db[which(db$abun%in%c(1,2,3,4,5,6) & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])
  th= glmSRnat.overall$impact.spread[i, "th.CI"]
  es <- as.numeric(glmSRnat.overall$est[i,])
  n <- as.numeric(glmSRnat.overall$n.obs[i,])[2:6]
  es[which(n<5)] <- NA
  es = which(!is.na(es))+1
  m <- max(es, na.rm=T)
  max.impact=as.character(db[which(db$abun>=6 & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])
  th.impact=as.character(db[which(db$abun>= th & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])

  plot(study_area, col="white", border = "grey60")
  plot(envplot[plot.alien, ],pch=22, cex = 0.5, col="tan",bg="tan", add=T)
  plot(envplot[th.impact, ],pch=22, cex = 0.5, col="sienna2", bg="sienna2", add=T)
  plot(envplot[max.impact, ],pch=22,  cex = 0.5, col="sienna4", bg="sienna4", add=T)

  mtext(3, text =sub("_", " ",species[i, "tip"]),adj = 0.1, line= -2, font=3, cex = 0.7)
}


legend(extent(study_area)@xmin,
       extent(study_area)@ymin , legend = c("< critical abundance",">critical abundance", "Dominant"),
       fill=c( "tan", "sienna2","sienna4"),xpd = NA, bty="n", border =c( "tan", "sienna2","sienna4"))


dev.off()


