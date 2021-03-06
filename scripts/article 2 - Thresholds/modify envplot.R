# modify envplot for GIS analysis:

coordinates(envplot) <- c("POINTX", "POINTY")
proj4string(envplot) <- original.CRS
envplot <- spTransform(envplot, CRS("+proj=longlat +datum=WGS84"))
envplot.grid <- points2grid(envplot, tolerance = 0.63, round=1)
grid <- SpatialGrid(envplot.grid, proj4string = proj4string(study_area))

# Extract landuse for point data in envplot
o = over(envplot,LUCAS)
envplot@data = cbind(envplot@data,o)
envplot$landuse <- envplot$PREV_LUC_N

## realgrasslands with landuse = lucasgrasslands  ( different from previous landcover data)
## dominant species
envplot$first.rank <-NA
firstranksp <- databp[databp$DominanceRank == 1 , c("SpeciesCode", "PlotName")]
firstranksp$ALIEN <- 0
firstranksp$ALIEN[firstranksp$SpeciesCode %in% aliens] <- 1
envplot$first.rank <-firstranksp$SpeciesCode[match(envplot$PLOTID,firstranksp$PlotName)]

herbaceouslands <- rownames(envplot@data)[which(
  envplot$first.rank %in% rownames(species)[species$Growth.forms %in% c('GR','HR')] )]

lucasgrasslands <- rownames(envplot@data)[which(
  (as.character(envplot$landuse) %in% c("Grassland - High producing",
                                        "Grassland - Low producing","Grassland - With woody biomass")))]

unimproved_lucasgrasslands <- rownames(envplot@data)[which(
  (as.character(envplot$landuse) %in% c("Grassland - Low producing","Grassland - With woody biomass")))]

realgrasslands <- lucasgrasslands[lucasgrasslands %in% herbaceouslands]
unimprovedgrasslands <- unimproved_lucasgrasslands[unimproved_lucasgrasslands %in% herbaceouslands]

#count occurrences of species
tmp <- colSums(comm[which(rownames(comm) %in% as.character(lucasgrasslands)),]>0, na.rm=T)
species$grassland.occur <- tmp[match(rownames(species), names(tmp))]

tmp <- colSums(comm[which(rownames(comm) %in% as.character(unimprovedgrasslands)),]>0, na.rm=T)
species$unimproved_grassland.occur <- tmp[match(rownames(species), names(tmp))]

## Transform ASPECT into northern orientation (cosinus) and Eatern (sinus):

envplot@data$Northern <- cospi(envplot@data$ASPECT/180)
envplot@data$NWestern <- cospi(envplot@data$ASPECT/180 - (1/4))
envplot@data$Eastern <- sinpi(envplot@data$ASPECT/180)

databp$Northern <-  cospi(databp$ASPECT/180)
databp$NWestern <- cospi(databp$ASPECT/180 - (1/4))
databp$Eastern <- sinpi(databp$ASPECT/180)
