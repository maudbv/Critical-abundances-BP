
### Mapping impacts

#############       Using GIS data
require(rgdal)
require(raster)
require(sp)
require(gstat)
library(raster)

study_area <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_study_area")
geology <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_geology_map")
population <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_population_data")
plots <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_points")

coordinates(envplot) <- c("POINTX", "POINTY")
proj4string(envplot) <- proj4string(study_area)

envplot.grid <- points2grid(envplot, tolerance = 0.63, round=1)
grid <- SpatialGrid(envplot.grid, proj4string = proj4string(study_area))


# as raster ?
envplot.ras <- raster(envplot)

### represent study area
### mapping landcover
landcover.col=c("grey","grey",
                "grey","grey",
                "grey","grey",
                "yellowgreen","grey",
                "forestgreen","grey",
                "grey","grey",
                "grey","grey",
                "grey","grey",
                "grey", "grey")
names(landcover.col) <- levels(as.factor(envplot$landcover))
landcover.col[c(2,4,8,11,12,13,14)] = "grey30"

plot(study_area, border = "grey", col="grey")
plot(envplot, add = T, pch=22, col = landcover.col[envplot$landcover],bg = landcover.col[envplot$landcover] )

legend('topright',legend= c("High productivity grassland", 
                            "Low productivity grassland","forest","other"),
       fill=c("goldenrod","forestgreen","red","grey"), cex=0.6)

#outliers in the NMDS.complete :
points(envplot[c("180", "787"),])


## proportion of invasion AR/TR

envplot$prop.alien <- envplot$SRali/envplot$SR

x <- cut(envplot$prop.alien , 20)
cp <- colorRampPalette(c("white" ,"tan","sienna4"))(20)
plot(study_area, border = "grey", col="grey")
plot(envplot, add = T, pch=22, col = cp[as.numeric(x)], bg = cp[as.numeric(x)])


### Impact spread :
pdf(file = "map impasct.pdf")
db=databp[databp$PlotName %in% realgrasslands,]
targets <- rownames(glmSRnat.overall$impact.spread)[
  which(!is.na(glmSRnat.overall$impact.spread$th.CI) &
          rownames(glmSRnat.overall$impact.spread)%in% aliens)]

par(mfrow=c(3,4), mar=c(1,1,1,1))                         
for( i in targets) {
  plot.alien=as.character(db[which(db$abun%in%c(1,2,3,4,5,6) & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])  
  th= glmSRnat.overall$impact.spread[i, "th.CI"]
  es <- as.numeric(glmSRnat.overall$est[i,])
  n <- as.numeric(glmSRnat.overall$n.obs[i,])[2:6]
  es[which(n<5)] <- NA
  es = which(!is.na(es))+1
  m <- max(es, na.rm=T)
  max.impact=as.character(db[which(db$abun>= m & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])
  th.impact=as.character(db[which(db$abun>= th & db$vegtype=="G" & db$SpeciesCode==i),"PlotName" ])

  plot(study_area, col="lightgrey", border = "lightgrey")
  plot(envplot[plot.alien, ],pch=22, cex = 0.5, col="tan",bg="tan", add=T)
  plot(envplot[th.impact, ],pch=22, cex = 0.5, col="sienna2", bg="sienna2", add=T)
  plot(envplot[max.impact, ],pch=22,  cex = 0.5, col="sienna4", bg="sienna4", add=T)
  plot(study_area, add =T, border= "grey60")
  mtext(3, text =sub("_", " ",species[i, "tip"]), line= -1, font=3, cex = 0.7)
}

plot.new()
legend("center", legend = c("< critical abundance",">critical abundance", "maximum abundance"),
       fill=c( "tan", "sienna2","sienna4"), bty="n", border =c( "tan", "sienna2","sienna4"))


dev.off()





### To be investigated for INTERPOLATION
library(sp)
library(lattice) # required for trellis.par.set():
trellis.par.set(sp.theme()) # sets color ramp to bpy.colors()

speciesab= envplot
sp = "TRIREP"
speciesab$abun <- db[db$SpeciesCode == sp,"abun"] [match(as.character(envplot$PLOTID),as.character(db[db$SpeciesCode == sp,"PlotName"]))] 
m

speciesab$abun [is.na(speciesab$abun)] =0

vSR.ok = variogram(log(abun)~1,speciesab)
ok.model = fit.variogram(vSR.ok, vgm(1, "Exp", 500, 1))
plot(vSR.ok, ok.model, main = "ordinary kriging")



abun.ok = krige(log(abun)~1, speciesab, grid, model = ok.model)
zn = abun.ok
zn[["a"]] = abun.ok[["var1.pred"]]
spplot(zn, c("a"), 
       names.attr = c("ordinary kriging"), 
       as.table = TRUE, main = "abundance interpolation")
       
       

alphaChannelSupported = function() { 
  !is.na(match(names(dev.cur()), c("pdf")))
}

data(meuse)
coordinates(meuse)=~x+y
data(meuse.riv)
meuse.sr = SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)),"meuse.riv")))
rv = list("sp.polygons", meuse.sr, 
          fill = ifelse(alphaChannelSupported(), "blue", "transparent"),
          alpha = ifelse(alphaChannelSupported(), 0.1, 1))
pts = list("sp.points", meuse, pch = 3, col = "grey", 
           alpha = ifelse(alphaChannelSupported(), .5, 1))
text1 = list("sp.text", c(180500,329900), "0", cex = .5, which = 4)
text2 = list("sp.text", c(181000,329900), "500 m", cex = .5, which = 4)
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(180500,329800), scale = 500, fill=c("transparent","black"), which = 4)

library(gstat, pos = match(paste("package", "sp", sep=":"), search()) + 1)
data(meuse.grid)
coordinates(meuse.grid) = ~x+y
gridded(meuse.grid) = TRUE
v.ok = variogram(log(zinc)~1, meuse)
ok.model = fit.variogram(v.ok, vgm(1, "Exp", 500, 1))
# plot(v.ok, ok.model, main = "ordinary kriging")
v.uk = variogram(log(zinc)~sqrt(dist), meuse)
uk.model = fit.variogram(v.uk, vgm(1, "Exp", 300, 1))
# plot(v.uk, uk.model, main = "universal kriging")
meuse[["ff"]] = factor(meuse[["ffreq"]])
meuse.grid[["ff"]] = factor(meuse.grid[["ffreq"]])
v.sk = variogram(log(zinc)~ff, meuse)
sk.model = fit.variogram(v.sk, vgm(1, "Exp", 300, 1))
# plot(v.sk, sk.model, main = "stratified kriging")
zn.ok = krige(log(zinc)~1,          meuse, meuse.grid, model = ok.model)
zn.uk = krige(log(zinc)~sqrt(dist), meuse, meuse.grid, model = uk.model)
zn.sk = krige(log(zinc)~ff,         meuse, meuse.grid, model = sk.model)
zn.id = krige(log(zinc)~1,          meuse, meuse.grid)

zn = zn.ok
zn[["a"]] = zn.ok[["var1.pred"]]
zn[["b"]] = zn.uk[["var1.pred"]]
zn[["c"]] = zn.sk[["var1.pred"]]
zn[["d"]] = zn.id[["var1.pred"]]
spplot(zn, c("a", "b", "c", "d"), 
       names.attr = c("ordinary kriging", "universal kriging with dist to river", 
                      "stratified kriging with flood freq", "inverse distance"), 
       as.table = TRUE, main = "log-zinc interpolation",
       sp.layout = list(rv, scale, text1, text2)
)

