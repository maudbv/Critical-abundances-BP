## Threshold master file

setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R_alienimpactBP/")
# setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R_alienimpactBP")


# load packages
library(doBy)
library(vegan)
library(markdown)
require(rgdal)
require(rgeos)
require(raster)
require(sp)
require(gstat)
require(maps)
library(raster)

### restore saved results:
# load("saved Rdata/article 2 - threshold/article threshold 1.2.Rdata")  ## all lucas grasslands
# load("saved Rdata/article 2 - threshold/article threshold 1.3.Rdata")  ## unimproved only lucas grasslands



### Import and modify data from scratch:
source('script/data/import BP species and environment data.R', encoding = "native.enc")
# source('script/data/import trait data.R', encoding = "native.enc")

# taxonomy solving :
# extract names, match them to DB, create a reference list of changes, update initial DB

# import GIS data
source('script/data/import GIS data.R', echo=FALSE)

# modify envplot for GIS analysis:

### add information from the envplot dataframe
coordinates(envplot) <- c("POINTX", "POINTY")
proj4string(envplot) <- proj4string(study_area)

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
unimprovedgrasslands <- lucasgrasslands[unimproved_lucasgrasslands %in% herbaceouslands]

tmp <- colSums(comm[which(rownames(comm) %in% as.character(lucasgrasslands)),]>0, na.rm=T)
species$grassland.occur <- tmp[match(rownames(species), names(tmp))]

tmp <- colSums(comm[which(rownames(comm) %in% as.character(unimprovedgrasslands)),]>0, na.rm=T)
species$unimproved_grassland.occur <- tmp[match(rownames(species), names(tmp))]

# load functions
source('script/functions/p2star.R')
source('script/article 2 - Thresholds/glm test.R')
source('script/article 2 - Thresholds/plotting functions.R')
source('script/article 2 - Thresholds/overall bootstrapped GLM.R')
source('script/article 2 - Thresholds/impact.size.R')
source('script/article 2 - Thresholds/summary.glm.R')
source('script/article 2 - Thresholds/correct th.CI.R')
source('script/article 2 - Thresholds/myraupcrick.R')
source('script/article 2 - Thresholds/dissim.nm.R')
source('script/functions/SMsim.R')
source('script/article 2 - Thresholds/Beta turnover and nestedness.R')
source('script/article 2 - Thresholds/gamma.trend.R')
source('script/article 2 - Thresholds/div partitioning.R')

### threshold analysis using GLMs

# set bootstrap sample size
nreps <- 999
 # db<- databp[databp$PlotName %in% realgrasslands,]
 db<- databp[databp$PlotName %in% unimprovedgrasslands,]
min.occur <- 5
min.class <- 2

#####  Bootstrpping within each class :

#  system.time(glmSR <- glm.test(db = db,bootstrap = T, min.occur= min.occur, nreps=nreps, CI=0.95, drastic =F,min.occur =min.occur,  min.class = min.occur))
#  glmSRnat <- glm.test(db = db,var="SRnat",bootstrap = T,min.occur= min.occur,  nreps=nreps, CI=0.95, drastic =F, min.occur =min.occur,  min.class = min.occur)
#  glmSRali <- glm.test(db = db,var="SRali",bootstrap = T,min.occur= min.occur,  nreps=nreps, CI=0.95, drastic =F,min.occur =min.occur,  min.class = min.occur)

# save(glmSR,glmSRnat,glmSRali, file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")
#load(file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")

#####  Overall bootstrapping :

# # create bootstrapped samples
#  system.time(boot.output <- bootstrap.dataset(db=db, min.occur =min.occur,  min.class = min.class, nreps = nreps))
#  system.time(boot.indices <- extract.indices(boot.output, db = db))
#
#  save( boot.output, boot.indices ,
#   file = "saved Rdata/article 2 - threshold/boot.output.2.2.Rdata")

## version 2.0 with 751 "realgrasslands"
## version 2.1 with 827 "lucasgrasslands" including the grasslands with woody biomass but still first ranked is herbaceous
## version 2.2 with 793 "unimprovedgrasslands" including the grasslands with woody biomass but still first ranked is herbaceous

# load(file = "saved Rdata/article 2 - threshold/boot.output.2.0.Rdata") # robin's landcover grasslands
 load(file = "saved Rdata/article 2 - threshold/boot.output.2.1.Rdata") # all lucas grasslands
# load(file = "saved Rdata/article 2 - threshold/boot.output.2.2.Rdata") # unimproved grasslands only

# calculate glms on bootstraps
# system.time(glmSR.overall <- glm.overallboot(db = db,boot.indices=boot.indices, sp.target = sp.target, variable = 'SR', min.occur= min.occur, min.class = min.class, nreps=nreps))
#  system.time(glmSRnat.overall <- glm.overallboot(db = db,boot.indices=boot.indices, sp.target = sp.target, variable = 'SRnat', min.occur= min.occur, min.class = min.class, nreps=nreps))
#  system.time(glmSRali.overall <- glm.overallboot(db = db,boot.indices=boot.indices, sp.target = sp.target, variable = 'SRali', min.occur= min.occur, min.class = min.class, nreps=nreps))
#
# save( glmSRnat.overall, glmSRali.overall,file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.2.Rdata")

# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.1.Rdata") # robin's landcover grasslands
 load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.1.Rdata") # all lucas grasslands
# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.2.Rdata") # unimproved grasslands only


### correct th.CI when miscalculated
glmSR.overall <- impact.spread(M = glmSR.overall, variable= "SR")
glmSRnat.overall <- impact.spread(M = glmSRnat.overall, variable= "SRnat")
glmSRali.overall <- impact.spread(M = glmSRali.overall, variable= "SRali")

### SUMMARY of frequencies of significantly negative effects and thresholds
glmSR.sum <-  summary.glmtest(M=glmSR.overall, group="ALIEN", type="overall.boot")
glmSRnat.sum <- summary.glmtest(M=glmSRnat.overall, group="ALIEN", type="overall.boot")
glmSRali.sum <- summary.glmtest(M=glmSRali.overall, group="ALIEN", type="overall.boot")


glmSR.sum.th <-  summary.glmtest(M=glmSR.overall, group="ALIEN", type="overall.boot",threshold= "th")
glmSRnat.sum.th <- summary.glmtest(M=glmSRnat.overall, group="ALIEN", type="overall.boot",threshold= "th")
glmSRali.sum.th <- summary.glmtest(M=glmSRali.overall, group="ALIEN", type="overall.boot", threshold= "th")

## overall bootstrap impact size
impact.SR <- impact.size (glmSR.overall)
impact.SRnat <- impact.size (glmSRnat.overall)
impact.SRali <- impact.size (glmSRali.overall)

 write.csv(impact.SRnat, "impact.SRnat.csv" )
 write.csv(impact.SRali, "impact.SRali.csv" )

 # Select pecies which show a threshold for Native richness
impsp <- rownames(glmSRnat.overall$impact.spread[which(!is.na(glmSRnat.overall$impact.spread$th.CI)
                                                        & (rownames(glmSRnat.overall$impact.spread) %in% aliens)),])

## alpha and beta trends
#  alpha.trend.nat <- alpha.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.rare", nreps = 999)
#  beta.trend.nat <- beta.trend(spnames = rownames(glmSRnat.overall$mean),null.model="permute.rare", nreps = 999)

###### gamma trends

# permute rare
 gamma.trend.nat <- gamma.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.rare", nreps = 999)

# permute all
 gamma.trend.nat.permute.all <- gamma.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.all", nreps = 999)

# resample beta all
gamma.trend.nat.beta <- gamma.trend(spnames = impsp, bootstrapped =TRUE, breps = 999,
                                                null.model="resample.beta", nreps = 999)
# resample beta rare
gamma.trend.nat.betarare <- gamma.trend(spnames = impsp, bootstrapped =TRUE, breps = 999,
                                    null.model="resample.beta.inrare", nreps = 999)

# # resample beta rare.vs.k  ### NOT WORKING
#  gamma.trend.nat.beta.rarevsk <- gamma.trend(spnames = impsp, bootstrapped =TRUE, breps = 99,
#                                      null.model="resample.beta.inrarevsk", nreps = 99)
#
#  save(gamma.trend.nat, alpha.trend.nat,gamma.trend.nat.permute.all,gamma.trend.nat.beta,gamma.trend.nat.betarare,
#  file = "saved Rdata/article 2 - threshold/gamma.trends.2.Rdata")
 save(gamma.trend.nat, alpha.trend.nat,gamma.trend.nat.permute.all,gamma.trend.nat.beta,gamma.trend.nat.betarare,
 file = "saved Rdata/article 2 - threshold/gamma.trends.unimproved.Rdata")

# load(file = "saved Rdata/article 2 - threshold/gamma.trends.Rdata")
 load(file = "saved Rdata/article 2 - threshold/gamma.trends.2.Rdata")
# load(file = "saved Rdata/article 2 - threshold/gamma.trends.unimproved.Rdata")

# ##### Gamma above/below #######
  # system.time(alpha.above.trend <- div.part.alpha.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "permute", nreps =999))

# # # permute total
#  system.time(gamma.above.nat.total <- div.part.gamma.nm(spnames = impsp, group=natives, null.model = "permute.total", nreps =99))
# #
# # # permute all
#   system.time(gamma.above.nat.permute.all <- div.part.gamma.nm(spnames = impsp, group=natives, null.model = "permute.all", nreps =999))
# #
# # # # resample all
#   system.time(gamma.above.nat.beta <- div.part.gamma.nm(spnames = impsp, group=natives, null.model = "resample.beta", nreps =99))
# #
# # # resample in rare
#   system.time(gamma.above.nat.betarare <- div.part.gamma.nm(spnames = impsp, group=natives, null.model = "resample.beta.inrare", nreps =999))
# #


# save(gamma.above.nat,gamma.above.nat.beta,gamma.above.nat.betarare,alpha.above.trend, file= "saved Rdata/article 2 - threshold/diversity.partitioning.Rdata")
#  save(alpha.above.trend,gamma.above.nat.total, gamma.above.nat.permute.all, gamma.above.nat.permute.all,gamma.above.nat.beta,gamma.above.nat.betarare,alpha.above.trend,       file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.Rdata")

 load (file= "saved Rdata/article 2 - threshold/diversity.partitioning.2.Rdata") # all lucas grasslands
# load (file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.Rdata") ##unimproved grasslands only


#### Beta diversity = turnover + nestedness   => NOT USED
 library(vegan)
# betasor.nat.all = betasor(spnames = "overall", group = natives, nreps = 999)
# betasor.ali.all = betasor(spnames = "overall", group = aliens, nreps = 999)
 betasor.nat = betasor(spnames = impsp, group = natives, nreps = 999)
 betasor.ali = betasor(spnames = impsp, group = aliens, nreps = 999)

 save(betasor.ali, betasor.nat, betasor.ali.all,betasor.nat.all,file="saved Rdata/article 2 - threshold/betasor results.Rdata")
 load(file="saved Rdata/article 2 - threshold/betasor results.Rdata")

#
# betasor.output <- as.data.frame(t(sapply(impsp, FUN = function(x) {
#   betasor_results <- t(betasor.nat[[x]][c("nestedness", "turnover", "sorensen"),
#                                         c("below.obs","above.obs","below.pval", "above.pval")])})))
# colnames(betasor.output) <- paste(rep(c("below.obs","above.obs","below.pval", "above.pval"),3),
#                                   c(rep("nestedness",4), rep("turnover",4), rep("sorensen",4)), sep = "_")
# betasor.output.ali <- as.data.frame(t(sapply(impsp, FUN = function(x) {
#   betasor_results <- t(betasor.ali[[x]][c("nestedness", "turnover", "sorensen"),
#                                         c("below.obs","above.obs","below.pval", "above.pval")])})))
# colnames(betasor.output.ali) <- paste(rep(c("below.obs","above.obs","below.pval", "above.pval"),3),
#                                   c(rep("nestedness",4), rep("turnover",4), rep("sorensen",4)), sep = "_")

## beta dissimilarities across grasslands

# community <- comm[which(rownames(comm) %in% realgrasslands),]
# community <- community[,colSums(community)>0]
# betadist.all <- betasor.dist(com = community)
#
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



# extract result table for loss in native diversity above thresholds:
sp.names = c(impsp)
out <- glmSRnat.overall$impact.spread
out <- out[sp.names,c("th.CI", "prevalence", "n.plot.impact")]
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

out$GRo <- sapply(sp.names, FUN = function(m) gamma.above.nat.total$gamma.below[m,out[m, "th"]] )
out$GRc <- sapply(sp.names, FUN = function(m) gamma.above.nat.total$gamma.above[m,out[m, "th"]] )
out$GRnull <- sapply(sp.names, FUN = function(m)  gamma.above.nat.total$null.above[m,out[m, "th"]]  )
out$GRsdnull <- sapply(sp.names, FUN = function(m)  gamma.above.nat.total$sdnull.above[m,out[m, "th"]]  )
out$gamma.loss <- out$GRnull - out$GRc

out$deltagamma <- sapply(sp.names, FUN = function(m)  gamma.above.nat$deltagamma[m,out[m, "th"]]  )
out$delta.z.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.nat$z.delta[m,out[m, "th"]]  )
out$delta.null.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.nat$null.delta[m,out[m, "th"]]  )
out$GRP.permute.all <- sapply(sp.names, FUN = function(m)  gamma.above.nat$P.delta[m,out[m, "th"]] )

out$delta.z.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$z.delta[m,out[m, "th"]]  )
out$delta.null.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$null.delta[m,out[m, "th"]]  )
out$GRP.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$P.delta[m,out[m, "th"]] )

out$BRo <- sapply(sp.names, FUN = function(m) gamma.above.nat$betap.below[m,out[m, "th"]] )
out$BRc <-sapply(sp.names, FUN = function(m) gamma.above.nat$betap.above[m,out[m, "th"]] )

# add betasor and beta nested => NOT RUN
 out=cbind(out, betasor.output[sp.names,])

# final formatting of table
table.div.part <- out[impsp,]

table.div.part <- orderBy(~ th.CI, table.div.part)

write.csv(data.frame(table.div.part), file="output diversity partitioning.csv")



# save results
 save.image("saved Rdata/article 2 - threshold/article threshold 1.2.Rdata")

 # save.image("saved Rdata/article 2 - threshold/article threshold 1.3.Rdata")

