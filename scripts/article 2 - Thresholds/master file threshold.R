## Threshold master file
setwd("/Users/maud/Documents/Work/Postdoc lincoln local /R/R_critical abundances analyses")

# load packages
library(doBy)
library(vegan)
library(markdown)
require(rgdal)
require(rgeos)
require(raster)
require(sp)
require(maps)
library(raster)
library(prettymapr)

### restore saved results: 
# load("saved Rdata/article 2 - threshold/article threshold 1.3.4.Rdata") 

# load functions   #############
source('scripts/functions/p2star.R')
source('scripts/article 2 - Thresholds/glm test.R')
source('scripts/article 2 - Thresholds/plotting functions.R')
source('scripts/article 2 - Thresholds/bootstrapping dataset.R')
# source('scripts/article 2 - Thresholds/overall bootstrapped GLM.R')
source('scripts/article 2 - Thresholds/overall bootstrapped GLM_with elevation.R')
source('scripts/article 2 - Thresholds/impact.size.R')
source('scripts/article 2 - Thresholds/summary.glm.R')
source('scripts/article 2 - Thresholds/correct th.CI.R')
source('scripts/article 2 - Thresholds/myraupcrick.R')
source('scripts/article 2 - Thresholds/dissim.nm.R')
source('scripts/functions/SMsim.R')
source('scripts/article 2 - Thresholds/Beta turnover and nestedness.R')
source('scripts/article 2 - Thresholds/gamma.trend.R')
source('scripts/article 2 - Thresholds/div partitioning.R')

### Import and modify data from scratch:  ############
source('scripts/data/import BP species and environment data.R', encoding = "native.enc")
# source('script/data/import trait data.R', encoding = "native.enc")

# import GIS data
source('scripts/data/import GIS data.R', echo=FALSE)

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

##### threshold analysis using GLMs  ########

# set bootstrap sample size
nreps <- 999
  # db<- databp[databp$PlotName %in% realgrasslands,]
db<- databp[databp$PlotName %in% unimprovedgrasslands,]
min.occur <- 5
min.class <- 2

#####  Bootstrpping within each class :

#  system.time(glmSR <- glm.test(db = db,bootstrap = T, nreps=nreps, CI=0.95, drastic =F,min.occur =min.occur,  min.class = min.occur))
#  glmSRnat <- glm.test(db = db,var="SRnat",bootstrap = T, nreps=nreps, CI=0.95, drastic =F, min.occur =min.occur,  min.class = min.occur)
  # glmSRali <- glm.test(db = db,var="SRali",bootstrap = T, nreps=nreps, CI=0.95, drastic =F,min.occur =min.occur,  min.class = min.occur)

# save(glmSR,glmSRnat,glmSRali, file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")
load(file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")

#####  Overall bootstrapping :

# # create bootstrapped samples
   # system.time(boot.output <- bootstrap.dataset(db=db, min.occur =min.occur,  min.class = min.class, nreps = nreps))
   # system.time(boot.indices <- extract.indices(boot.output, db = db))
# #
   # save( boot.output, boot.indices ,
   #  file = "saved Rdata/article 2 - threshold/boot.output.2.3.Rdata")

## version 2.0 with 751 "realgrasslands"
## version 2.1 with 827 "lucasgrasslands" including the grasslands with woody biomass but still first ranked is herbaceous
## version 2.2 with 595 "unimprovedgrasslands" including the grasslands with woody biomass but still first ranked is herbaceous

# load(file = "saved Rdata/article 2 - threshold/boot.output.2.0.Rdata") # robin's landcover grasslands
# load(file = "saved Rdata/article 2 - threshold/boot.output.2.1.Rdata") # all lucas grasslands
# load(file = "saved Rdata/article 2 - threshold/boot.output.2.2.Rdata") # unimproved grasslands only
 load(file = "saved Rdata/article 2 - threshold/boot.output.2.3.Rdata") # unimproved grasslands only - RUN 2
   
# calculate glms on bootstraps
# system.time(glmSR.overall <- glm.overallboot(db = db,boot.indices=boot.indices, sp.target = sp.target, variable = 'SR', min.occur= min.occur, min.class = min.class, nreps=nreps))
 # system.time(glmSRnat.overall <- glm.overallboot(db = db,boot.ind =boot.indices, variable = 'SRnat',covar = c("DEM_10","SLOPE", "Northern", "SRali"), min.occur= min.occur, min.class = min.class, nreps=nreps))
 # system.time(glmSRali.overall <- glm.overallboot(db = db,boot.indices=boot.indices, sp.target = sp.target, variable = 'SRali', min.occur= min.occur, min.class = min.class, nreps=nreps))
#
 # save( glmSRnat.overall,file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.5.with elevation+northern+srali.Rdata")

# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.1.Rdata") # robin's landcover grasslands
# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.1.Rdata") # all lucas grasslands
# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.2.Rdata") # unimproved grasslands only
   
   load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.2.Rdata") ## Old version without cofctors
   glmSRnat.overall.nocovar <- glmSRnat.overall
   
# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.4.with elevation+northern+srali.Rdata") # unimproved grasslands only + elevation as cofactor + "covar.tab" name forgotten in output
 load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.5.with elevation+northern+srali.Rdata") # unimproved grasslands only + elevation as cofactor + bootindices RUN 2
      
### Modify impact.spread calculations when necessary
# glmSR.overall <- correct.impact.spread(M = glmSR.overall, variable= "SR")
glmSRnat.overall <- correct.impact.spread(M = glmSRnat.overall, db=db, variable= "SRnat",threshold.type= "") # strictly negative and significant indices following a critical abundance
glmSRnat.overall <- correct.impact.spread(M = glmSRnat.overall, db=db, variable= "SRnat",threshold.type= "custom") 

# glmSRali.overall <- correct.impact.spread(M = glmSRali.overall,db=db, variable= "SRali")

### SUMMARY of frequencies of significantly negative effects and thresholds
# glmSR.sum <-  summary.glmtest(M=glmSR.overall, group="ALIEN", type="overall.boot")
glmSRnat.sum <- summary.glmtest(M=glmSRnat.overall,   group="ALIEN", type="overall.boot")
# glmSRali.sum <- summary.glmtest(M=glmSRali.overall, group="ALIEN", type="overall.boot")


# glmSR.sum.th <-  summary.glmtest(M=glmSR.overall, db=db,group="ALIEN", type="overall.boot",threshold= "th")
glmSRnat.sum.th <- summary.glmtest(M=glmSRnat.overall,group="ALIEN", type="overall.boot",threshold= "th")
# glmSRali.sum.th <- summary.glmtest(M=glmSRali.overall, group="ALIEN", type="overall.boot", threshold= "th")

## overall bootstrap impact size
# impact.SR <- impact.size (glmSR.overall)
impact.SRnat <- impact.size (glmSRnat.overall)
# impact.SRali <- impact.size (glmSRali.overall)

 write.csv(impact.SRnat, "impact.SRnat.csv" )
 # write.csv(impact.SRali, "impact.SRali.csv" )

 
 
 # Select pecies which show a threshold for Native richness
impsp <- rownames(glmSRnat.overall$impact.spread[which(!is.na(glmSRnat.overall$impact.spread$th.CI) & (rownames(glmSRnat.overall$impact.spread) %in% aliens)),])
# 
# impsp.ali <- rownames(glmSRali.overall$impact.spread[which(!is.na(glmSRali.overall$impact.spread$th.CI)),])

## alpha and beta trends
#  alpha.trend.nat <- alpha.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.rare", nreps = 999)
#  beta.trend.nat <- beta.trend(spnames = rownames(glmSRnat.overall$mean),null.model="permute.rare", nreps = 9)


##### gamma trends ###############
# # permute rare
# gamma.trend.nat <- gamma.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.rare", nreps = 999)
#
# # permute all
# gamma.trend.nat.permute.all <- gamma.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.all", nreps = 999)

# # resample beta all
# gamma.trend.nat.beta <- gamma.trend(spnames = impsp, bootstrapped =TRUE, breps = 999,
                                                # null.model="resample.beta", nreps = 999)
 # # resample beta rare
# gamma.trend.nat.betarare <- gamma.trend(spnames = impsp, bootstrapped =TRUE, breps = 999,
#                                     null.model="resample.beta.inrare", nreps = 999)

# # resample beta rare.vs.k  ### NOT WORKING
#  gamma.trend.nat.beta.rarevsk <- gamma.trend(spnames = impsp, bootstrapped =TRUE, breps = 99,
#                                      null.model="resample.beta.inrarevsk", nreps = 99)
#
#  save(gamma.trend.nat, alpha.trend.nat,gamma.trend.nat.permute.all,gamma.trend.nat.beta,gamma.trend.nat.betarare,
#  file = "saved Rdata/article 2 - threshold/gamma.trends.2.Rdata")
#   save(gamma.trend.nat, alpha.trend.nat,gamma.trend.nat.permute.all,gamma.trend.nat.beta,
#   file = "saved Rdata/article 2 - threshold/gamma.trends.unimproved.Rdata")

# save(gamma.trend.nat,   file = "saved Rdata/article 2 - threshold/gamma.trends.unimproved.2.1.Rdata")

# load(file = "saved Rdata/article 2 - threshold/gamma.trends.Rdata")
# load(file = "saved Rdata/article 2 - threshold/gamma.trends.2.Rdata")
load(file = "saved Rdata/article 2 - threshold/gamma.trends.unimproved.2.1.Rdata")


##### Gamma above/below #######
# system.time(alpha.above.trend <- div.part.alpha.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "permute", nreps =999))
#
## permute total
# system.time(gamma.above.nat.total <- div.part.gamma.nm(spnames =rownames(glmSRnat.overall$mean), group=natives, null.model = "permute.total", nreps =99))
#
# permute all
# system.time(gamma.above.nat.permute.all <- div.part.gamma.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "permute.all", nreps =999))
#
## resample all
# system.time(gamma.above.nat.beta <- div.part.gamma.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "resample.beta", nreps =999))
#
## resample in rare
# system.time(gamma.above.nat.betarare <- div.part.gamma.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "resample.beta.inrare", nreps =999))
#
#  save(alpha.above.trend, gamma.above.nat.permute.all,
#     gamma.above.nat.beta,gamma.above.nat.betarare,
#    file= "saved Rdata/article 2 - threshold/diversity.partitioning.Rdata")
#
#  save(alpha.above.trend,
#    gamma.above.nat.permute.all,gamma.above.nat.beta,
#    file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.Rdata")
#
  # save(alpha.above.trend, gamma.above.nat.permute.all,
  #   file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.withcofactors.Rdata")
 
# load (file= "saved Rdata/article 2 - threshold/diversity.partitioning.Rdata") # all lucas grasslands
 load (file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.withcofactors.Rdata") ##unimproved grasslands only



##### gamma trends ALIEN richness ###############

#  # permute rare
#   gamma.trend.ali <- gamma.trend(spnames = rownames(glmSRali.overall$mean), group = aliens,
#                                  null.model="permute.rare", nreps = 999)
#   save(gamma.trend.ali,
#   file = "saved Rdata/article 2 - threshold/gamma.trends.SRali.Rdata")
#  load(file = "saved Rdata/article 2 - threshold/gamma.trends.SRali.Rdata")
#
#

# ##### Gamma and alpha above/below  ALIEN richness #######
#
#  system.time(alpha.above.trend.ali <- div.part.alpha.nm(spnames = rownames(glmSRali.overall$mean), group=aliens, null.model = "permute", nreps =999))
# #permute all
# system.time(gamma.above.ali.permute.all <- div.part.gamma.nm(spnames =  rownames(glmSRali.overall$mean), group=aliens, null.model = "permute.all", nreps =999))
#
# #resample all
# system.time(gamma.above.ali.beta <- div.part.gamma.nm(spnames = rownames(glmSRali.overall$mean), group=aliens, null.model = "resample.beta", nreps =99))
#
# save(alpha.above.trend.ali, gamma.above.ali.permute.all,
# gamma.above.ali.beta, file= "saved Rdata/article 2 - threshold/diversity.partitioning.SRali.Rdata")
#
# load( file= "saved Rdata/article 2 - threshold/diversity.partitioning.SRali.Rdata")



##### Beta diversity = turnover + nestedness   => USED   ###########
#  library(vegan)
    # betasor.nat = betasor.multi(spnames = impsp, group = natives, null.comm = T, bootstrap = F, nreps = 999)
    # betasor.boot.nat = betasor.multi(spnames = impsp, group = natives, null.comm = F, bootstrap = T, nreps = 999)
    # 

# #  save(betasor.ali, betasor.nat, betasor.ali.all,betasor.nat.all,file="saved Rdata/article 2 - threshold/betasor results.Rdata")
# #  load(file="saved Rdata/article 2 - threshold/betasor results.Rdata")
# 
# save(betasor.nat,betasor.boot.nat, file = "saved Rdata/beta diversity output.2.1.Rdata") #with cofactors and RUN2

 load(file = "saved Rdata/beta diversity output.2.1.Rdata")
  
table.sorensen <-  as.data.frame(t(data.frame(lapply(betasor.nat, FUN = function(x) x[3,c("obs.below", "obs.above", "obs.diff", "z.diff", "p.diff")]))))
table.sorensen.boot <-  as.data.frame(t(data.frame(lapply(betasor.nat.boot, FUN = function(x) x[3,c("obs.below", "mean.below","sd.below", "p.below",
                                                                                                    "obs.above","mean.above","sd.above", "p.above",
                                                                                                    "obs.diff", "mean.diff","sd.diff", "p.diff")]))))


table.turnover <- t(data.frame(lapply(betasor.nat, FUN = function(x) x[1,c("obs.below", "obs.above", "obs.diff", "z.diff", "p.diff")])))
table.nestedness <- as.data.frame(t(data.frame(lapply(betasor.nat, FUN = function(x) x[2,c("obs.below", "obs.above", "obs.diff", "z.diff", "p.diff")]))))
table.nestedness.boot <-  as.data.frame(t(data.frame(lapply(betasor.nat.boot, FUN = function(x) x[2,c("obs.below", "mean.below","sd.below", "p.below",
                                                                                                      "obs.above","mean.above","sd.above", "p.above",
                                                                                                      "obs.diff", "mean.diff","sd.diff", "p.diff")]))))


## beta dissimilarities across grasslands

# community <- comm[which(rownames(comm) %in% realgrasslands),]
# community <- community[,colSums(community)>0]
# betadist.all <- betasor.dist(com = community)

# community <- comm[which(rownames(comm) %in% realgrasslands),]
# community <- community[,colSums(community)>0]
# community<- community[,colnames(community) %in% natives]
# betadist.nat <- betasor.dist(com = community)
#
# community <- comm[which(rownames(comm) %in% realgrasslands),]
# community <- community[,colSums(community)>0]
# community<- community[,colnames(community) %in% aliens]
# betadist.ali <- betasor.dist(com = community)

#### SPATIAL CLUSTERING
# source('C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R - alienimpactBP/script/article 2 - Thresholds/spatial correlation.R')
# save(spDst, spatial.mnnd, mnnd.th, mnnd.trends, mnnd.above,
#      file = "saved Rdata/article 2 - threshold/spatial MNND.Rdata"))

# load(file = "saved Rdata/article 2 - threshold/spatial MNND.Rdata")


  
  

# extract result table for loss in native diversity above thresholds:   ################
sp.names = impsp
out <- glmSRnat.overall$impact.spread
out <- out[sp.names,c("th.CI", "prevalence", "n.plot.impact", "n.plot.dominant")]
out<- cbind(species = species[sp.names, "tip"], out)

out$th.CI.SRali <- glmSRali.overall$impact.spread[sp.names,"th.CI"]
out$th <- out$th.CI
out$th[is.na(out$th)] <- out$th.CI.SRali[is.na(out$th)]

out$aRc.old <- sapply(sp.names, FUN = function(m) glmSRnat.overall$mean[m,out[m, "th"]] )
out$aRo.old <- glmSRnat.overall$mean[sp.names, "C1"]

out$aRo <- sapply(sp.names, FUN = function(m) alpha.above.trend$alpha.below[m,out[m, "th"]] )
out$aRc <- sapply(sp.names, FUN = function(m) alpha.above.trend$alpha.above[m,out[m, "th"]] )
out$aRo.sd <- sapply(sp.names, FUN = function(m) alpha.above.trend$alpha.below.sd[m,out[m, "th"]] )
out$aRc.sd <- sapply(sp.names, FUN = function(m) alpha.above.trend$alpha.above.sd[m,out[m, "th"]] )

out$aRnull <- sapply(sp.names, FUN = function(m) alpha.above.trend$null.above[m,out[m, "th"]]  )
out$aRsdnull <- sapply(sp.names, FUN = function(m)  alpha.above.trend$sdnull.above[m,out[m, "th"]]  )
out$alpha.loss <- out$aRnull  - out$aRc
out$aR.P <- sapply(sp.names, FUN = function(m)  alpha.above.trend$P.above[m,out[m, "th"]] )

out$deltaalpha <- sapply(sp.names, FUN = function( m) alpha.above.trend$deltaalpha[m,out[m, "th"]]  )
out$deltaalpha.null <- sapply(sp.names, FUN = function(m)  alpha.above.trend$null.delta[m,out[m, "th"]]  )
out$deltaalpha.z.permute.all <- sapply(sp.names, FUN = function(m) {
  (alpha.above.trend$deltaalpha[m,out[m, "th"]] - alpha.above.trend$null.delta[m,out[m, "th"]])/alpha.above.trend$sdnull.delta[m,out[m, "th"]]
  })

out$aR.t <- sapply(sp.names, FUN = function(m)  alpha.above.trend$t[m,out[m, "th"]] )
out$aR.df<- sapply(sp.names, FUN = function(m)  alpha.above.trend$df[m,out[m, "th"]] )
out$aR.Pt <- sapply(sp.names, FUN = function(m)  alpha.above.trend$P.t[m,out[m, "th"]] )

out$GRo <- sapply(sp.names, FUN = function(m) gamma.above.nat.permute.all$gamma.below[m,out[m, "th"]] )
out$GRc <- sapply(sp.names, FUN = function(m) gamma.above.nat.permute.all$gamma.above[m,out[m, "th"]] )
out$GRnull <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$null.above[m,out[m, "th"]]  )
out$GRsdnull <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$sdnull.above[m,out[m, "th"]]  )
out$gamma.loss <- out$GRnull - out$GRc

out$deltagamma <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$deltagamma[m,out[m, "th"]]  )
out$delta.z.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$z.delta[m,out[m, "th"]]  )
out$delta.null.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$null.delta[m,out[m, "th"]]  )
out$GRP.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.nat.permute.all$P.delta[m,out[m, "th"]] )

out$deltagamma.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$deltagamma[m,out[m, "th"]]  )
out$delta.z.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$z.delta[m,out[m, "th"]]  )
out$delta.null.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$null.delta[m,out[m, "th"]]  )
out$GRP.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$P.delta[m,out[m, "th"]] )

out$BRo <- sapply(sp.names , FUN = function(m) gamma.above.nat.permute.all$betap.below[m,out[m, "th"]] )
out$BRc <-sapply(sp.names, FUN = function(m) gamma.above.nat.permute.all$betap.above[m,out[m, "th"]] )

out$z.beta.diff <-table.turnover[sp.names,"z.diff"]
out$P.beta.diff <-table.turnover[sp.names,"p.diff"]

# final formatting of table
table.div.part <- out[impsp,]


write.csv(as.matrix(orderBy(~ th.CI, table.div.part)), file="output diversity partitioning.csv", fileEncoding= "native.enc")


# ## Table 2 :  divpart for alien richness

# # extract result table for loss in alien diversity above thresholds:   ################
# sp.names = impsp.ali
# out.ali <- glmSRali.overall$impact.spread
# out.ali <- out.ali[sp.names,c("th.CI", "prevalence", "n.plot.impact")]
# out.ali<- cbind(species = species[sp.names, "tip"], out.ali)
#
# out.ali$th.CI.SRali <- glmSRali.overall$impact.spread[sp.names,"th.CI"]
# out.ali$th <- out.ali$th.CI
# out.ali$th[is.na(out.ali$th)] <- out.ali$th.CI.SRali[is.na(out.ali$th)]
#
# out.ali$aRc.old <- sapply(sp.names, FUN = function(m) glmSRali.overall$mean[m,out.ali[m, "th"]] )
# out.ali$aRo.old <- glmSRali.overall$mean[sp.names, "C1"]
#
# out.ali$aRo <- sapply(sp.names, FUN = function(m) alpha.above.trend.ali $alpha.below[m,out.ali[m, "th"]] )
# out.ali$aRc <- sapply(sp.names, FUN = function(m) alpha.above.trend.ali $alpha.above[m,out.ali[m, "th"]] )
# out.ali$aRo.sd <- sapply(sp.names, FUN = function(m) alpha.above.trend.ali $alpha.below.sd[m,out.ali[m, "th"]] )
# out.ali$aRc.sd <- sapply(sp.names, FUN = function(m) alpha.above.trend.ali $alpha.above.sd[m,out.ali[m, "th"]] )
#
#
# out.ali$aRnull <- sapply(sp.names, FUN = function(m) alpha.above.trend.ali $null.above[m,out.ali[m, "th"]]  )
# out.ali$aRsdnull <- sapply(sp.names, FUN = function(m)  alpha.above.trend.ali $sdnull.above[m,out.ali[m, "th"]]  )
# out.ali$alpha.loss <- out.ali$aRnull  - out.ali$aRc
# out.ali$aR.P <- sapply(sp.names, FUN = function(m)  alpha.above.trend.ali $P.above[m,out.ali[m, "th"]] )
#
# out.ali$deltaalpha <- sapply(sp.names, FUN = function( m) alpha.above.trend.ali $deltaalpha[m,out.ali[m, "th"]]  )
# out.ali$deltaalpha.null <- sapply(sp.names, FUN = function(m)  alpha.above.trend.ali $null.delta[m,out.ali[m, "th"]]  )
# out.ali$deltaalpha.z.permute.all <- sapply(sp.names, FUN = function(m) {
#   (alpha.above.trend.ali $deltaalpha[m,out.ali[m, "th"]] - alpha.above.trend.ali $null.delta[m,out.ali[m, "th"]])/alpha.above.trend.ali $sdnull.delta[m,out.ali[m, "th"]]
# })
#
# out.ali$aR.t <- sapply(sp.names, FUN = function(m)  alpha.above.trend.ali $t[m,out.ali[m, "th"]] )
# out.ali$aR.df<- sapply(sp.names, FUN = function(m)  alpha.above.trend.ali $df[m,out.ali[m, "th"]] )
# out.ali$aR.Pt <- sapply(sp.names, FUN = function(m)  alpha.above.trend.ali $P.t[m,out.ali[m, "th"]] )
#
# out.ali$GRo <- sapply(sp.names, FUN = function(m) gamma.above.ali.permute.all$gamma.below[m,out.ali[m, "th"]] )
# out.ali$GRc <- sapply(sp.names, FUN = function(m) gamma.above.ali.permute.all$gamma.above[m,out.ali[m, "th"]] )
# out.ali$GRnull <- sapply(sp.names, FUN = function(m)  gamma.above.ali.permute.all$null.above[m,out.ali[m, "th"]]  )
# out.ali$GRsdnull <- sapply(sp.names, FUN = function(m)  gamma.above.ali.permute.all$sdnull.above[m,out.ali[m, "th"]]  )
# out.ali$gamma.loss <- out.ali$GRnull - out.ali$GRc
#
# out.ali$deltagamma <- sapply(sp.names, FUN = function(m)  gamma.above.ali.permute.all$deltagamma[m,out.ali[m, "th"]]  )
# out.ali$delta.z.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.ali.permute.all$z.delta[m,out.ali[m, "th"]]  )
# out.ali$delta.null.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.ali.permute.all$null.delta[m,out.ali[m, "th"]]  )
# out.ali$GRP.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.ali.permute.all$P.delta[m,out.ali[m, "th"]] )
#
# out.ali$deltagamma.beta <- sapply(sp.names, FUN = function(m)  gamma.above.ali.beta$deltagamma[m,out.ali[m, "th"]]  )
# out.ali$delta.z.beta <- sapply(sp.names, FUN = function(m)  gamma.above.ali.beta$z.delta[m,out.ali[m, "th"]]  )
# out.ali$delta.null.beta <- sapply(sp.names, FUN = function(m)  gamma.above.ali.beta$null.delta[m,out.ali[m, "th"]]  )
# out.ali$GRP.beta <- sapply(sp.names, FUN = function(m)  gamma.above.ali.beta$P.delta[m,out.ali[m, "th"]] )
#
# # final formatting of table
# table.div.part.ali <- out.ali[impsp.ali,]
#
# table.div.part.ali <- orderBy(~ th.CI, table.div.part.ali)
#
# write.csv(as.matrix(table.div.part.ali), file="output diversity partitioning.SRali.csv", fileEncoding= "native.enc")






# ### Table S1 with AR and NR critical abundances  ##############
# sp.names= rownames(glmSRnat.overall$impact.spread)[!is.na(glmSRnat.overall$impact.spread$th.CI)  |
#                                                      !is.na(glmSRali.overall$impact.spread$th.CI)]
# 
# out <- glmSRnat.overall$impact.spread
# out <- out[sp.names,c("prevalence", "th.CI",  "n.plot.impact")]
# out.ali <- glmSRali.overall$impact.spread
# out.ali <- out.ali[sp.names,c("th.CI", "n.plot.impact")]
# 
# out<- cbind(species = species[sp.names, c("SpeciesName", "ALIEN")], out, out.ali)
# 
# names(out) <- c("SpeciesName","ALIEN", "prevalence", "th.nat", "nplot.nat", "th.ali", "nplot.ali")
# 
# out$aRc.nat.old <- sapply(sp.names, FUN = function(m) {
#   if (!is.na(out[m, "th.nat"])) glmSRnat.overall$mean[m,out[m, "th.nat"]] else NA
#   })
# out$aRo.nat.old <- glmSRnat.overall$mean[sp.names, "C1"]
# 
# out$aRc.ali.old <- sapply(sp.names, FUN = function(m) {
#   if (!is.na(out[m, "th.ali"])) glmSRali.overall$mean[m,out[m, "th.ali"]] else NA
# })
# out$aRo.ali.old <- glmSRnat.overall$mean[sp.names, "C1"]
# 
# write.csv(as.matrix(out), file="critical abundances AR and NR.csv", fileEncoding= "native.enc")


####### TABLE 1 GLMS (only significant species) :
table1 <- data.frame( nobs = glmSRnat.overall$impact.spread[impsp,]$prevalence,
                      glmSRnat.overall$glms[impsp,],
                      DEM10.c = glmSRnat.overall$covar.tab$DEM_10[impsp,]$coef.glm,
                      DEM10.P = glmSRnat.overall$covar.tab$DEM_10[impsp,]$P.coef,
                      
                      SLOPE.c = glmSRnat.overall$covar.tab$SLOPE[impsp,]$coef.glm,
                      SLOPE.P = glmSRnat.overall$covar.tab$SLOPE[impsp,]$P.coef,

                      Northness.c = glmSRnat.overall$covar.tab$Northern[impsp,]$coef.glm,
                      Northness.P = glmSRnat.overall$covar.tab$Northern[impsp,]$P.coef,

                      SRali.c = glmSRnat.overall$covar.tab$SRali[impsp,]$coef.glm,
                      SRali.P = glmSRnat.overall$covar.tab$SRali[impsp,]$P.coef,
                    
                      abun.meanc = rowMeans(glmSRnat.overall$est[impsp,], na.rm = T)
)


table1 <- table1 [order(table.div.part$th.CI),]
write.csv(table1, file= "table1.csv", row.names = row.names(table1))

####### TABLE S2 GLMS (ALL focal species):

tableS2init <- data.frame( nobs = glmSRnat.overall.nocovar$impact.spread$prevalence,
                       glmSRnat.overall.nocovar$glms,
                       abun.meanc = rowMeans(glmSRnat.overall.nocovar$est, na.rm = T)
                       ,
                       abun.minP = apply(glmSRnat.overall.nocovar$P,1, min, na.rm = T)
)



tableS2 <- data.frame( nobs = glmSRnat.overall$impact.spread$prevalence,
                      glmSRnat.overall$glms,
                      DEM10.c = glmSRnat.overall$covar.tab$DEM_10$coef.glm,
                      DEM10.P = glmSRnat.overall$covar.tab$DEM_10$P.coef,
                      
                      SLOPE.c = glmSRnat.overall$covar.tab$SLOPE$coef.glm,
                      SLOPE.P = glmSRnat.overall$covar.tab$SLOPE$P.coef,
                      
                      Northness.c = glmSRnat.overall$covar.tab$Northern$coef.glm,
                      Northness.P = glmSRnat.overall$covar.tab$Northern$P.coef,
                      
                      SRali.c = glmSRnat.overall$covar.tab$SRali$coef.glm,
                      SRali.P = glmSRnat.overall$covar.tab$SRali$P.coef,
                      
                      abun.meanc = exp(rowMeans(glmSRnat.overall$est, na.rm = T))
)

write.csv(tableS2, file= "tableS2.csv", row.names = row.names(tableS1))


####### TABLE 2 :
table2 <- data.frame(Acrit = table.div.part$th.CI,
                     Presence = table.div.part$prevalence,
                     above.Acrit = table.div.part$n.plot.impact,
                     dominance = table.div.part$n.plot.dominant,
                     alpha.below= round(table.div.part$aRo, 2),
                     alpha.above = round(table.div.part$aRc,2),
                     delta.alpha = round(table.div.part$aRc -  table.div.part$aRo, 2),
                     perc.delta.alpha = round(((table.div.part$aRc -  table.div.part$aRo)/table.div.part$aRo)*100, 1),
                     gamma.below = table.div.part$GRo,
                     gamma.above = table.div.part$GRc,
                     delta.gamma.c = round(table.div.part$deltagamma - table.div.part$delta.null.permute.all,2),
                     delta.gamma.P = round(table.div.part$GRP.permute.all, 4),
                      delta.beta.z = table.turnover[rownames(table.div.part),"z.diff"],
                      delta.beta.P = table.turnover[rownames(table.div.part),"p.diff"]
                     
)
rownames(table2) <- table.div.part$species
table2 <- table2 [order(table.div.part$th.CI),]
write.csv(table2, file= "table2.csv", row.names = row.names(table2))





# save results   ############
# save.image("saved Rdata/article 2 - threshold/article threshold 1.3.4.Rdata")

