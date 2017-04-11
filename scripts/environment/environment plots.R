## Amalyse environmental conditions in data set

require(ade4)
require(FactoMineR)
library(raster)

# distribution des especes entre G et W ?

nativesinG=na.omit(unique(databp$SpeciesCode[databp$vegtype=="G" & databp$ALIEN==0]))
length(nativesinG)/365
nativesinW=na.omit(unique(databp$SpeciesCode[databp$vegtype=="W" & databp$ALIEN==0]))
length(nativesinW)/365
nativesinboth=na.omit(nativesinG[which(nativesinG %in% nativesinW)])
length(nativesinboth)/365
nativesG=na.omit(nativesinG[which(!nativesinG %in% nativesinW)])
length(nativesG)/365
nativesW=na.omit(nativesinW[which(!nativesinW %in% nativesinG)])
length(nativesW)/365
length(nativesinboth)+length(nativesG)+length(nativesW)

aliensinG=na.omit(unique(databp$SpeciesCode[databp$vegtype=="G" & databp$ALIEN==1]))
length(aliensinG)/365
length(aliensinG)/538
aliensinW=na.omit(unique(databp$SpeciesCode[databp$vegtype=="W" & databp$ALIEN==1]))
length(aliensinW)/487

aliensinboth=na.omit(aliensinG[which(aliensinG %in% aliensinW)])
length(aliensinboth)/308
aliensG=na.omit(aliensinG[which(!aliensinG %in% aliensinW)])
length(aliensG)/308
aliensW=na.omit(aliensinW[which(!aliensinW %in% aliensinG)])
length(aliensW)/308
length(aliensinboth)+length(aliensG)+length(aliensW)


length(which(envplot$SRnat==0 & envplot$SRali>0 ))  # 161
length(which(envplot$SRnat==1 & envplot$SRnat>0 ))  # 121
(161+121)/1124



# create a Raster MAP for the peninsula
makeraster=function(var=envplot$DEM_10, coordX="newX", coordY="newY", data=envplot){
  spg <- data.frame(var, data[, c(coordX,coordY )])
  names(spg)=c("var","x", "y")
  coordinates(spg) <- ~ x + y
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  return(rasterDF)
}

a=makeraster(var=as.numeric(!is.na(envplot$PLOTID)))
b=makeraster(var=envplot$DEM_10)
c=makeraster(var=as.numeric(envplot$vegtype=="W")+1)
d=makeraster(var=envplot$SRnat)
BPraster=stack(a,b,c,d)
names(BPraster)=c("PLOTID","DEM_10","vegtype", "SRnat")

x11()
par(mar=c(1,1,1,1))
plot(BPraster$DEM_10, col="goldenrod",legend=F, axes=F)
plot(BPraster$DEM_10, col=colorRampPalette(colors=c("grey80","grey10"))(225),alpha=0.8, add=T, legend=F, axes=F)
plot(BPraster$vegtype, col=colorRampPalette(c("goldenrod","forestgreen"))(2), add=T, alpha=0.4, legend=F, axes=F)


### Ordination of env parameters
            
tmp=envplot[,c(  "PANN0080","PSUM0080", "PWIN0080",  "AVTEMP0080","MINTEM0080",
                     "MAXTEM0080","GDD","SOLRADR","MOISTR","UR1991_DNS",
                     "BLDG_DIST","DIST_HWSLD","DIST_METAL","DIST_PRIV","DIST_SRIV" ,
                     "DEM_10","ASPECT","SLOPE","PH_MID","SR")]
tmp=tmp-colMeans(tmp)/ apply(tmp, 2, sd, na.rm=T)    
tmp=cbind(tmp,envplot[,c("SRnat","SRali", "vegtype")])

## PCA on restricted variables
bb=PCA(tmp[,c( "PANN0080", "AVTEMP0080","SOLRADR","MOISTR","DEM_10","PH_MID",
                     "ASPECT","SLOPE","SRnat","SRali", "vegtype")],
       quanti.sup=c(9,10), quali.sup=11)
# 
# bb=PCA(tmp[,c( "PANN0080", "AVTEMP0080","DEM_10",
#              "ASPECT","SLOPE","SRnat","SRali", "vegtype")],
#        quanti.sup=c(6,7), quali.sup=8)
# 
# bb=PCA(tmp[,c("DEM_10","ASPECT","SLOPE","SRnat","SRali", "vegtype")],
#        quanti.sup=c(4,5), quali.sup=6)


plot(bb, choix="ind",habillage=11,col.hab=c('goldenrod', 'forestgreen'),
     label="none", axes=c(1,2), ylim=c(-6,6), xlim=c(-2,2), xaxt='n', yaxt='n', ann=F)
# aa <- cbind.data.frame(envplot$vegtype,bb$ind$coord)
# ell <- coord.ellipse(aa,bary=TRUE,level.conf = 0.99)
points(bb$ind$coord[as.character(weirdwood),1:2])
par(new=T)
plot(bb, choix="ind", label="none",col.ind=NA, axes=c(1,2),
     col.quali='black', ylim=c(-6,6), xlim=c(-2,2), xaxt='n', yaxt='n', ann=F)
par(new=T)
plot(bb, choix="var", axes=c(1,2), ylim=c(-1,1), xlim=c(-2,2))


plot(bb, choix="ind",habillage=11,col.hab=c('goldenrod', 'forestgreen'),
     label="none", axes=c(1,3),  ylim=c(-6,6), xlim=c(-2,2), xaxt='n', yaxt='n', ann=F)
par(new=T)               
plot(bb, choix="var", axes=c(1,3), ylim=c(-1,1), xlim=c(-2,2))

#### vegtype effect on env ?
summary(lm(SLOPE ~ vegtype, data=tmp))
summary(lm(MOISTR ~ vegtype, data=tmp))
summary(lm(DEM_10 ~ vegtype, data=tmp))

#### GRASSLANDS environmental factors
envplot.grass=envplot[envplot$PLOTID%in% grasslands,]
bb=PCA(envplot.grass[,c( "PANN0080", "AVTEMP0080","SOLRADR","MOISTR","DEM_10",
                     "ASPECT","SLOPE","SRnat","SRali")],quanti.sup=c(8,9))
                        
plot(bb, choix="ind", label="none", col.ind="goldenrod",axes=c(1,2), ylim=c(-6,6), xlim=c(-2,2), xaxt='n', yaxt='n', ann=F)
par(new=T)
plot(bb, choix="var", axes=c(1,2), ylim=c(-1,1), xlim=c(-2,2))

plot(bb, choix="ind", label="none", col.ind="goldenrod",axes=c(1,3), ylim=c(-6,6), xlim=c(-2,2), xaxt='n', yaxt='n', ann=F)
par(new=T)
plot(bb, choix="var", axes=c(1,3), ylim=c(-1,1), xlim=c(-2,2))

     
########map of vegetation types
plot(envplot$POINTX,envplot$POINTY,pch=20, col=c('goldenrod', 'forestgreen')[as.numeric(envplot$vegtype=="W")+1])
points(envplot[as.character(weirdwood), "POINTX"],envplot[as.character(weirdwood), "POINTY"], pch=20, col="red")
######## Plots dominated by aliens             
library(geometry)
library(tripack)

envplot.grass=tmp[envplot$PLOTID%in% grasslands,]
bb=PCA(envplot.grass[,c( "PANN0080", "AVTEMP0080","SOLRADR","MOISTR","DEM_10",
                     "ASPECT","SLOPE","SRnat","SRali")],quanti.sup=c(8,9))

dist.envir=dist(bb$ind$coord[, 1:3])   


########## Environmental ordination + convex hull volume
par(mfrow=c(5,4), mar=c(1,2,3,1))                         
for( i in rownames(effect.summary)[which(effect.summary$total.effect>0 & rownames(effect.summary)%in% aliens)] )
  {
plot.alien=as.character(databp[which(databp$domlevels%in%c(1,2,3,4,5) & databp$vegtype=="G" & databp$SpeciesCode==i),"PlotName" ]   )

plot(bb, choix="ind", label="none", col.ind="goldenrod",axes=c(1,2),
    xaxt='n', yaxt='n', ann=F)

# H=convhulln(bb$ind$coord[plot.alien, c(1,2)])
H=chull(bb$ind$coord[plot.alien, c(1,2)])
polygon(bb$ind$coord[plot.alien, c(1,2)][H,], border="blue")

H=tri.mesh(x=bb$ind$coord[plot.alien,1],y=bb$ind$coord[plot.alien, 2])
inHull=which(in.convex.hull(H, x=bb$ind$coord[, 1], y=bb$ind$coord[, 2]))
points(bb$ind$coord[inHull, c(1,2)], col="blue", pch=20)
     title(main=i, line=1)
points(bb$ind$coord[plot.alien, c(1,2)], col="firebrick", pch=20)
     title(main=i, line=1)

}
                         
# plot(bb, choix="ind", label="none", col.ind="goldenrod",axes=c(1,3), ylim=c(-6,6), xlim=c(-2,2), xaxt='n', yaxt='n', ann=F)
# points(bb$ind$coord[plot.lolium, c(1,3)], col="firebrick", pch=20)
# par(new=T)
# plot(bb, choix="var", axes=c(1,3), ylim=c(-1,1), xlim=c(-2,2))
# 

########## convex hull volume on raw data

E=envplot.grass[,c( "PANN0080", "AVTEMP0080","SOLRADR","MOISTR","DEM_10",
                     "ASPECT","SLOPE","SRnat","SRali")]
E=apply(E,MARGIN=2,FUN=function(x) sapply(x, FUN=function(y) (mean(x, na.rm=T)-y)/sd(x)) )
i="ACHMIL"

par(mfrow=c(5,4), mar=c(1,2,3,1))                         
for( i in rownames(effect.summary)[which(effect.summary$total.effect>0 & rownames(effect.summary)%in% aliens)] )
  {
plot.alien=as.character(databp[which(databp$domlevels%in%c(1,2,3,4,5) & databp$vegtype=="G" & databp$SpeciesCode==i),"PlotName" ]   )

plot(bb, choix="ind", label="none", col.ind="goldenrod",axes=c(1,2),
    xaxt='n', yaxt='n', ann=F)

H=convhulln(E[1:100,1:3], options="p")

polygon(bb$ind$coord[plot.alien, c(1,2)][H,], border="blue")

H=tri.mesh(x=bb$ind$coord[plot.alien,1],y=bb$ind$coord[plot.alien, 2])
inHull=which(in.convex.hull(H, x=bb$ind$coord[, 1], y=bb$ind$coord[, 2]))
points(bb$ind$coord[inHull, c(1,2)], col="blue", pch=20)
     title(main=i, line=1)
points(bb$ind$coord[plot.alien, c(1,2)], col="firebrick", pch=20)
     title(main=i, line=1)

}


## Mapping of points
par(mfrow=c(5,4), mar=c(1,2,3,1))                         
for( i in rownames(effect.summary)[which(effect.summary$total.effect>0 & rownames(effect.summary)%in% aliens)] )
  {
plot.alien=as.character(databp[which(databp$domlevels%in%c(1,2,3,4,5) & databp$vegtype=="G" & databp$SpeciesCode==i),"PlotName" ]   )

plot(envplot$POINTX,envplot$POINTY,pch=20, col=c('goldenrod', 'white')[as.numeric(envplot$vegtype=="W")+1])
points(envplot[plot.alien, c("POINTX", "POINTY")],pch=20, col="firebrick")
}

# mapping landcover Robin


landcover.col=c("brown","forestgreen",
                "grey","green",
                "orange","grey",
                "goldenrod","olivedrab",
                "palegoldenrod","blue",
                "pink","purple",
                "black","darkgrey",
                "slateblue","yellow",
                "grey", "red")

landcover.col=c("grey","grey",
                "grey","grey",
                "red","grey",
                "goldenrod","grey",
                "forestgreen","grey",
                "grey","grey",
                "grey","grey",
                "grey","grey",
                "grey", "grey")
quartz()
par(mar=c(1,1,1,1))
 plot(envplot$POINTX,envplot$POINTY,pch=22,xlim=c(min(envplot$POINTX), 2550000),
        col=landcover.col[as.numeric(as.factor(envplot$landcover))], bg=landcover.col[as.numeric(as.factor(envplot$landcover))])
  
# legend('topright',legend= levels(as.factor(envplot$landcover)), fill=landcover.col, cex=0.6)
legend('topright',legend= c("High productivity grassland", 
                            "Low productivity grassland", "Gorse and Broom","other"),
       fill=c("goldenrod","forestgreen","red","grey"), cex=0.6)

[1] "Afforestation (imaged, post LCDB 1)" "Broadleaved Indigenous Hardwoods"   
[3] "Built-up Area"                       "Forest Harvested"                   
[5] "Gorse and Broom"                     "Grey Scrub"                         
[7] "High Producing Exotic Grassland"     "Indigenous Forest"                  
[9] "Low Producing Grassland"             "Major Shelterbelts"                 
[11] "Manuka and or Kanuka"                "Other Exotic Forest"                
[13] "Pine Forest - Closed Canopy"         "Pine Forest - Open Canopy"          
[15] "River and Lakeshore Gravel and Rock" "Short-rotation Cropland"            
[17] "Urban Parkland/ Open Space"          "Vineyard"   

### plot species i am looking for
par(mfrow=c(3,3), mar=c(1,2,3,1))                         
for( i in c("GNAAUD","JUNAUS","LINBIE", "RUMPUL","DICCRI","RYTUNA", "COTAUS", "CENTEN")) {
 plot.alien=as.character(databp[which(databp$abun%in%c(1,2,3,4,5,6) & databp$vegtype=="G" & databp$SpeciesCode==i),"PlotName" ]   )
plot(envplot$POINTX,envplot$POINTY,pch=20, 
     main=i,
     col=c('lightgrey', 'olivedrab')[as.numeric(envplot$vegtype=="W")+1])
points(envplot[plot.alien, c("POINTX", "POINTY")],pch=20, col="firebrick")
}


                         