<<<<<<< HEAD
### Mapping impacts


### mapping landcover
landcover.col=c("grey","grey",
                "grey","grey",
                "red","grey",
                "goldenrod","grey",
                "forestgreen","grey",
                "grey","grey",
                "grey","grey",
                "grey","grey",
                "grey", "grey")
par(mar=c(1,1,1,1))
plot(envplot$POINTX,envplot$POINTY,pch=22,xlim=c(min(envplot$POINTX), 2550000),
     col=landcover.col[as.numeric(as.factor(envplot$landcover))], bg=landcover.col[as.numeric(as.factor(envplot$landcover))])
legend('topright',legend= c("High productivity grassland", 
                            "Low productivity grassland", "Gorse and Broom","other"),
       fill=c("goldenrod","forestgreen","red","grey"), cex=0.6)


### plot species i am looking for
db=databp[databp$PlotName %in% realgrasslands,]
targets <- rownames(glmSRnat$boot.thresh)[which(!is.na(glmSRnat$boot.thresh$th))]

par(mfrow=c(4,4), mar=c(1,2,3,1))                         
for( i in targets) {
   plot.alien=as.character(db[which(db$abun%in%c(1,2,3,4,5,6) & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])  
  th= glmSRnat$boot.thresh[i, "th"]
  impact=as.character(db[which(db$abun>= th & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])
  
  plot(envplot$POINTX,envplot$POINTY,pch=22,main=i,col='lightgrey', bg="lightgrey",  axes=F)
  points(envplot[plot.alien, c("POINTX", "POINTY")],pch=22, col="goldenrod",bg="goldenrod")
  points(envplot[impact, c("POINTX", "POINTY")],pch=22,  col="firebrick", bg="firebrick")
}

=======
### Mapping impacts


### mapping landcover
landcover.col=c("grey","grey",
                "grey","grey",
                "red","grey",
                "goldenrod","grey",
                "forestgreen","grey",
                "grey","grey",
                "grey","grey",
                "grey","grey",
                "grey", "grey")
par(mar=c(1,1,1,1))
plot(envplot$POINTX,envplot$POINTY,pch=22,xlim=c(min(envplot$POINTX), 2550000),
     col=landcover.col[as.numeric(as.factor(envplot$landcover))], bg=landcover.col[as.numeric(as.factor(envplot$landcover))])
legend('topright',legend= c("High productivity grassland", 
                            "Low productivity grassland", "Gorse and Broom","other"),
       fill=c("goldenrod","forestgreen","red","grey"), cex=0.6)


### plot species i am looking for
db=databp[databp$PlotName %in% realgrasslands,]
targets <- rownames(glmSRnat.grass$boot.thresh)[which(!is.na(glmSRnat.grass$boot.thresh$th))]

par(mfrow=c(3,5), mar=c(1,2,3,1))                         
for( i in targets) {
   plot.alien=as.character(db[which(db$abun%in%c(1,2,3,4,5,6) & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])  
  th= glmSRnat.grass$boot.thresh[i, "th"]
  impact=as.character(db[which(db$abun>= th & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])
  
  plot(envplot$POINTX,envplot$POINTY,pch=22,main=i,col='lightgrey', bg="lightgrey",  axes=F)
  points(envplot[plot.alien, c("POINTX", "POINTY")],pch=22, col="olivedrab",bg="olivedrab")
  points(envplot[impact, c("POINTX", "POINTY")],pch=22,  col="firebrick", bg="firebrick")
}

>>>>>>> Post-doc-lincoln/master
