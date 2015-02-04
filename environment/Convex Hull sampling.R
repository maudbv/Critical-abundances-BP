# Volume Subsampling of plots and testing invaded vs. non invaded plots
library(geometry)
library(tripack)
library(hypervolume)
library(raster)
library(FactoMineR)

## Mapping of points
par(mfrow=c(5,4), mar=c(1,2,3,1))                         
for( i in rownames(effect.summary)[which(effect.summary$total.effect>0 & rownames(effect.summary)%in% aliens)] )
{
  plot.alien=as.character(databp[which(databp$domlevels%in%c(1,2,3,4,5) & databp$vegtype=="G" & databp$SpeciesCode==i),"PlotName" ]   )
  
  plot(envplot$newX,envplot$newY,pch=20, col=c('goldenrod', 'white')[as.numeric(envplot$vegtype=="W")+1])
  points(envplot[plot.alien, c("newX", "newY")],pch=20, col="firebrick")
}


########## HYPERVOLUME on raw data for grassland aliens

# standardized environemntal values for grassland plots
E=envplot.grass[,c( "PANN0080", "AVTEMP0080","SOLRADR","DEM_10","ASPECT","SLOPE")]
E=apply(E,MARGIN=2,FUN=function(x) sapply(x, FUN=function(y) (mean(x, na.rm=T)-y)/sd(x)) )

# Example of total hypervolume of all grassland plots
bw=estimate_bandwidth(E)
H=hypervolume(E,repsperpoint=1000, bandwidth=bw, quantile=0.1, name="HV_allplots")
plot(H) 

# Creating list of subsamples of plots for each "important" alien 
######! (Will have to be extended to include more aliens)
alien.names=rownames(effect.summary)[which(effect.summary$total.effect>0 & rownames(effect.summary)%in% aliens)]
HValienlist=list()
subsampleslist=list()

for( i in 1:length(alien.names)  )  {
#i="ACHMIL"
plot.alien=as.character(databp[which(databp$domlevels%in%c(1,2,3,4,5) & databp$vegtype=="G" & databp$SpeciesCode==alien.names[i]),"PlotName" ]   )
bw=estimate_bandwidth(E[plot.alien,])
HValienlist[[i]]=hypervolume(E[plot.alien,],repsperpoint=1000, bandwidth=bw, quantile=0.1, name=paste("HV",alien.names[i],sep="_"))
subsample=hypervolume_inclusion_test(HValienlist[[i]], E, reduction_factor = 1, verbose = T)
subsampleslist[[i]]=E[subsample,]
}
# save(HValienlist,subsampleslist,alien.names, file="saved Rdata/HValien.Rdata")
     
# Mapping the subsamples within hypervolumes
x11()
par(mfrow=c(5,4), mar=c(1,2,3,1)) 
for( i in 1:length(alien.names)  )  {  
# Mapping on banks peninsula
plot.alien=as.character(databp[which(databp$domlevels%in%c(1,2,3,4,5) & databp$vegtype=="G" & databp$SpeciesCode==alien.names[i]),"PlotName" ]   ) 
plot(envplot$POINTX,envplot$POINTY,pch=22, axes=F,cex=0.3,
     col=c('goldenrod', 'grey')[as.numeric(envplot$vegtype=="W")+1],
     bg=c('goldenrod', 'grey')[as.numeric(envplot$vegtype=="W")+1])
points(envplot[rownames(subsampleslist[[i]]), c("POINTX", "POINTY")],cex=0.3,pch=22, col="forestgreen",bg="forestgreen")
points(envplot[plot.alien, c("POINTX", "POINTY")],pch=22,cex=0.3, col="firebrick",bg="firebrick")
title(main=alien.names[i])
}

# Ordination of the subsamples within hypervolumes

envplot.grass=envplot[envplot$PLOTID%in% grasslands,]
bb=PCA(envplot.grass[,c( "PANN0080", "AVTEMP0080","SOLRADR","DEM_10",
                         "ASPECT","SLOPE","SRnat","SRali")],quanti.sup=c(8,9))
plot3d(bb$ind$coord[,1:3])
dist.envir=dist(bb$ind$coord[, 1:3])   

x11()
par(mfrow=c(5,4), mar=c(1,2,3,1)) 
for( i in 1:length(alien.names)  )  {  
  # Mapping on banks peninsula
  plot.alien=as.character(databp[which(databp$domlevels%in%c(1,2,3,4,5) & databp$vegtype=="G" & databp$SpeciesCode==alien.names[i]),"PlotName" ]   ) 
  plot(bb, choix="ind", label="none", col.ind="grey",axes=c(1,2),
       xaxt='n', yaxt='n', ann=F)
  
  points(bb$ind$coord[rownames(subsampleslist[[i]]),1:2],pch=20, col="orange")
  points(bb$ind$coord[plot.alien,1:2],pch=20, col="red")
  
  title(main=alien.names[i])
}



########## Environmental ordination + convex hull volume
i=10
plot.alien=as.character(databp[which(databp$domlevels%in%c(1,2,3,4,5) & databp$vegtype=="G" & databp$SpeciesCode==alien.names[i]),"PlotName" ]   )
   
  tmp=bb$ind$coord[plot.alien, c(1:3)]
  
  plot3d(bb$ind$coord[, c(1:3)],  col="goldenrod", xaxt='n', yaxt='n', ann=F) 
  H=t(convhulln(tmp))
  rgl.triangles(x=tmp[H,1], y=tmp[H,2],z=tmp[H,3], col = "#0000FF", add=T, alpha=0.2)

#   H=tri.mesh(x=bb$ind$coord[plot.alien,1],y=bb$ind$coord[plot.alien, 2])
#   inHull=which(in.convex.hull(H, x=bb$ind$coord[, 1], y=bb$ind$coord[, 2]))
#   points(bb$ind$coord[inHull, c(1,2)], col="blue", pch=20)
#   title(main=i, line=1)
#   points(bb$ind$coord[plot.alien, c(1,2)], col="firebrick", pch=20)
#   title(main=i, line=1)
  
