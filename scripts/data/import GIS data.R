
#############       Import GIS data and update envplot object with land use data
require(rgdal)
require(rgeos)
require(raster)
require(sp)
require(maps)
library(raster)
require (maptools)

# import Federico data
study_area <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_study_area")
plots <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_points")
# geology <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_geology_map")
# population <- readOGR(dsn = "data/GIS", layer ="Banks_Peninsula_population_data")

# Homogenize CRS
original.CRS <- proj4string(study_area)
study_area <-  spTransform(study_area, CRS("+proj=longlat +datum=WGS84"))
plots <-  spTransform(plots, CRS("+proj=longlat +datum=WGS84"))
# population <-  spTransform(population, CRS("+proj=longlat +datum=WGS84"))
# geology <-  spTransform(geology, CRS("+proj=longlat +datum=WGS84"))


### import landuse classification of new zealand (LUCAS)
LUCAS <- readOGR(dsn = "data/GIS", layer ="lucas-bp")
LUCAS <-  spTransform(LUCAS, CRS("+proj=longlat +datum=WGS84"))
LUCAS <- crop(LUCAS, extent(study_area) +  c(-0.02, +0.02, -.02, +.02))
region <- unionSpatialPolygons(LUCAS, IDs = rep(1,length(LUCAS)), threshold=1)


## Map of study area
plot(region, col = "grey50", border = NA)
polygon(extent(region)[c(1,2,2,1)], extent(region)[c(3,3,4,4)], col ="lightblue", border = NA)
plot(region, add= T, col = "grey70", border = NA)
plot(study_area, add= T, col = "grey90", border = NA)
plot(plots, add= T,  col ="goldenrod",bg ="goldenrod", pch = 22, cex = 0.7)
detach("package:GISTools", unload=TRUE)
library(maps)
map.scale(ratio = FALSE)


# ## colors for land uses
# LUCAS.col =c( "violet","violetred2","chartreuse4","olivedrab1","goldenrod","sienna4","antiquewhite3",
#               "steelblue4","black", "white","skyblue")
# # [1] "Cropland - Annual"              "Cropland - Perennial"           "Grassland - High producing"
# # [4] "Grassland - Low producing"      "Grassland - With woody biomass" "Natural Forest"
# # [7] "Other"                          "Planted Forest - Pre-1990"      "Settlements"
# # [10] "Wetland - Open water"           "Wetland - Vegetated non forest"
# names(LUCAS.col) = levels(LUCAS@data$PREV_LUC_N)
# 
# 
# # Map of land use and survey
# par(mar=c(0,0,2,2))
# plot(LUCAS, border=LUCAS.col[LUCAS@data$PREV_LUC_N],col= LUCAS.col[LUCAS@data$PREV_LUC_N])
# points(envplot[realgrasslands,], pch =3, col="black")
# legend("topright", legend=names(LUCAS.col), fill = LUCAS.col, cex=0.7)
# 
# scalebar <- SpatialPoints(rbind(p1 <-c(bbox(LUCAS)[1,1],bbox(LUCAS)[2,1]),
#                                 p2 <- c(bbox(LUCAS)[1,1] + 5000,bbox(LUCAS)[2,1])),
#                           proj4string =CRS(proj4string(study_area)))
# 
# segments(bbox(LUCAS)[1,1],bbox(LUCAS)[2,1], bbox(LUCAS)[1,1] + 5000,bbox(LUCAS)[2,1], lwd=1.2)
# text(  bbox(LUCAS)[1,1] + 2500,bbox(LUCAS)[2,1]+1000,label= "5 km")
# 
# 
# tmp <- list(x=2523922, y=5717542)
# arrows(tmp$x,tmp$y, tmp$x,tmp$y+2500, code=2, length=0.1 , lwd=1.2)
# text(  tmp$x,tmp$y+3250,label= "N")


