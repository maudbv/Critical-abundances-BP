## Threshold master file
setwd("/Users/maud/Documents/Work/Postdoc lincoln local /R/R_critical abundances analyses")

# load packages 
library(doBy)
library(vegan)
library(rgdal)
library(rgeos)
library(raster)
library(sp)
library(maps)
library(raster)
library(data.table)


# load functions   #############
source('scripts/functions/p2star.R')
source('scripts/article 2 - Thresholds/glm test.R')
source('scripts/article 2 - Thresholds/plotting functions.R')
source('scripts/article 2 - Thresholds/bootstrapping dataset.R')
# source('scripts/article 2 - Thresholds/overall bootstrapped GLM.R')
# source('scripts/article 2 - Thresholds/overall bootstrapped GLM_with elevation.R')
source('scripts/article 2 - Thresholds/overall bootstrapped GLM_with covariables.R')
source('scripts/article 2 - Thresholds/impact.size.R')
source('scripts/article 2 - Thresholds/summary.glm.R')
source('scripts/article 2 - Thresholds/correct th.CI.with positives.R')
source('scripts/article 2 - Thresholds/myraupcrick.R')
source('scripts/article 2 - Thresholds/dissim.nm.R')
source('scripts/functions/SMsim.R')
source('scripts/article 2 - Thresholds/Beta turnover and nestedness.R')
source('scripts/article 2 - Thresholds/gamma.trend.R')
source('scripts/article 2 - Thresholds/div partitioning.R')


# restore saved results:  ####
# load("saved Rdata/article 2 - threshold/article threshold 1.3.7.Rdata") 

#### Import and modify data from scratch:  ############
# source('scripts/data/import BP species and environment data.R', encoding = "native.enc")
load('saved Rdata/Banks_Peninsula_data.Rdata')

# source('script/data/import trait data.R', encoding = "native.enc")

# import GIS data
# source('scripts/data/import GIS data.R', echo=FALSE)
load('saved Rdata/GIS data.Rdata')

# update envplot
source('scripts/article 2 - Thresholds/modify envplot.R')

#save.image(file = "saved Rdata/Critical_abundance_base_data.Rdata")

#### Threshold analysis using GLMs  ########

# set bootstrap sample size
nreps <- 999
  # db<- databp[databp$PlotName %in% realgrasslands,]
db<- databp[databp$PlotName %in% unimprovedgrasslands,]
min.occur <- 5
min.class <- 2

#####  Bootstrpping within each class (OLD) :________________________________________________________________

#  glmSRnat <- glm.test(db = db,var="SRnat",bootstrap = T, nreps=nreps, CI=0.95, drastic =F, min.occur =min.occur,  min.class = min.occur)
  # glmSRali <- glm.test(db = db,var="SRali",bootstrap = T, nreps=nreps, CI=0.95, drastic =F,min.occur =min.occur,  min.class = min.occur)

# save(glmSR,glmSRnat,glmSRali, file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")
# load(file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")

#####  Overall bootstrapping :_____________________________________________________________________________

# # create bootstrapped samples
# system.time(boot.output <- bootstrap.dataset(db=db, min.occur =min.occur,  min.class = min.class, nreps = nreps))
# system.time(boot.indices <- extract.indices(boot.output, db = db))
# save( boot.output, boot.indices,file = "saved Rdata/article 2 - threshold/boot.output.2.3.Rdata")

## version 2.0 with 751 "realgrasslands"
## version 2.1 with 827 "lucasgrasslands" including the grasslands with woody biomass but still first ranked is herbaceous
## version 2.2 with 595 "unimprovedgrasslands" including the grasslands with woody biomass but still first ranked is herbaceous

# load(file = "saved Rdata/article 2 - threshold/boot.output.2.0.Rdata") # robin's landcover grasslands
# load(file = "saved Rdata/article 2 - threshold/boot.output.2.1.Rdata") # all lucas grasslands
 load(file = "saved Rdata/article 2 - threshold/boot.output.2.2.Rdata") # unimproved grasslands only
# load(file = "saved Rdata/article 2 - threshold/boot.output.2.3.Rdata") # unimproved grasslands only - RUN 2
   
# ________calculate glms on bootstraps _____________________________________________________________________
 
# NO covariables:
# system.time(glmSRnat.overall <- glm.overallboot(db = db,boot.ind =boot.indices, variable = 'SRnat',
#                                                 covar = c("DEM_10","SLOPE", "Northern", "SRali"),
#                                                 min.occur= min.occur, min.class = min.class, nreps=nreps))
# 
# system.time(glmSRali.overall <- glm.overallboot(db = db,boot.indices=boot.indices, variable = 'SRali',
#                                                 covar = c("DEM_10","SLOPE", "Northern", "SRali"),
#                                                 min.occur= min.occur, min.class = min.class, nreps=nreps))
# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.1.Rdata") # robin's landcover grasslands
# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.1.Rdata") # all lucas grasslands
# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.2.Rdata") # unimproved grasslands only
load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.2.Rdata") ## Old version without cofctors
glmSRnat.overall.nocovar <- glmSRnat.overall  # Store a version with no covariables
rm (glmSRnat.overall)

# Covariables for the JoAE resubmission (OLD):
# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.5.with elevation+northern+srali.Rdata") 
 # unimproved grasslands only + elevation as cofactor + bootindices RUN 2
      
# # Covariables + distance to bldg *************TO BE RE-RUN ************:
# system.time(glmSRnat.bldg <- glm.overallboot(db = db,boot.ind =boot.indices, variable = 'SRnat',
#                                                 covar = c("DEM_10","SLOPE", "Northern","BLDG_DIST", "SRali"),
#                                                 min.occur= min.occur, min.class = min.class, nreps=nreps))
# system.time(glmSRali.bldg <- glm.overallboot(db = db,boot.ind =boot.indices, variable = 'SRnat',
#                                                 covar = c("DEM_10","SLOPE", "Northern","BLDG_DIST", "SRali"),
#                                                 min.occur= min.occur, min.class = min.class, nreps=nreps))
 # save( glmSRnat.bldg,
 #      file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.8_with distance to roads.Rdata")
  load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.8_with distance to roads.Rdata")

  
## RESUBMISSION to JoE April 2018
## ALL Covariables included in GLMS (without distance to building)  *************TO BE RUN ************:
# system.time(glmSRnat<- glm.overallboot(db = db,boot.ind =boot.indices, variable = 'SRnat',
#                                                 covar = c("DEM_10","SLOPE", "Northern","SRali"),
#                                                 min.occur= min.occur, min.class = min.class, nreps=nreps))
# system.time(glmSRali <- glm.overallboot(db = db,boot.ind =boot.indices, variable = 'SRali',
#                                                 covar = c("DEM_10","SLOPE", "Northern"),
#                                                 min.occur= min.occur, min.class = min.class, nreps=nreps))
# save( glmSRnat, glmSRali,
#       file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.9_all covariables except bldg.Rdata")
load( file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.9_all covariables except bldg.Rdata")
 glmSRnat.covar <- glmSRnat
 glmSRali.covar <- glmSRali
 
## Major Revisions for JoE April 2018 
## ALL Covariables included in GLMS (without distance to building)  *************TO BE RUN ************:
 # system.time(glmSRnat<- glm.overallboot(db = db,boot.ind =boot.indices, variable = 'SRnat',
 #                                                covar = c("year","DEM_10","SLOPE", "Northern","SRali"),
 #                                                 min.occur= min.occur, min.class = min.class, nreps=nreps))
# 
#  system.time(glmSRali<- glm.overallboot(db = db,boot.ind =boot.indices, variable = 'SRali',
#                                                 covar = c("year","DEM_10","SLOPE", "Northern","SRali"),
#                                                  min.occur= min.occur, min.class = min.class, nreps=nreps))
# 
#  
 # save( glmSRnat, glmSRali,
 #       file = "saved Rdata/article 2 - threshold/overall.boot.glms.3_all covariables include year.Rdata")

 load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.3_all covariables include year.Rdata")


### Calculate impact.spread and critical abundances (= "threshold type")  ####

glmSRnat.overall.withbldg <- correct.impact.spread(M = glmSRnat.bldg, db=db, variable= "SRnat",threshold.type= "custom")
glmSRnat.overall.nocovar <- correct.impact.spread(M = glmSRnat.overall.nocovar, db=db, variable= "SRnat",threshold.type= "custom") 
glmSRnat.covar <- correct.impact.spread(M = glmSRnat.covar, db=db, variable= "SRnat",threshold.type= "custom") 


glmSRnat.overall<- correct.impact.spread(M = glmSRnat, db=db, variable= "SRnat",threshold.type= "custom") 
glmSRali.overall<- correct.impact.spread(M = glmSRali, db=db, variable= "SRali",threshold.type= "custom") 

### SUMMARY of frequencies of significantly negative effects and thresholds
glmSRnat.sum <- summary.glmtest(M=glmSRnat.overall,   group="ALIEN", type="overall.boot")
glmSRali.sum <- summary.glmtest(M=glmSRali.overall, group="ALIEN", type="overall.boot")

glmSRnat.covar.sum <- summary.glmtest(M=glmSRnat.covar,   group="ALIEN", type="overall.boot")
glmSRnat.nocovar.sum <- summary.glmtest(M=glmSRnat.overall.nocovar,   group="ALIEN", type="overall.boot")
glmSRnat.bldg.sum <- summary.glmtest(M=glmSRnat.overall.withbldg,   group="ALIEN", type="overall.boot")


## overall bootstrap impact size
impact.SRnat <- impact.size (glmSRnat.overall)
impact.SRali <- impact.size (glmSRali.overall)

impact.SRnat.nocovar <- impact.size (glmSRnat.overall.nocovar)


 write.csv(impact.SRnat, "impact.SRnat.csv" )
 write.csv(impact.SRali, "impact.SRali.csv" )

 
 
 # Select pecies which show a threshold for Native richness
impsp <- rownames(glmSRnat.overall$impact.spread[which(!is.na(glmSRnat.overall$impact.spread$th.CI) & (rownames(glmSRnat.overall$impact.spread) %in% aliens)),])

impsp.nocovar <- rownames(glmSRnat.overall.nocovar$impact.spread[which(!is.na(glmSRnat.overall.nocovar$impact.spread$th.CI) & (rownames(glmSRnat.overall.nocovar$impact.spread) %in% aliens)),])
 
impsp.ali <- rownames(glmSRali.overall$impact.spread[which(!is.na(glmSRali.overall$impact.spread$th.CI)),])

# important positives

impsp.pos <- rownames(glmSRnat.overall$impact.spread[which(!is.na(glmSRnat.overall$impact.spread$pth.CI) & (rownames(glmSRnat.overall$impact.spread) %in% aliens)),])


impsp.pos.nocovar <- rownames(glmSRnat.overall.nocovar$impact.spread[which(!is.na(glmSRnat.overall.nocovar$impact.spread$pth.CI) & (rownames(glmSRnat.overall.nocovar$impact.spread) %in% aliens)),])

####  alpha and beta trends
#  alpha.trend.nat <- alpha.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.rare", nreps = 999)
#  beta.trend.nat <- beta.trend(spnames = rownames(glmSRnat.overall$mean),null.model="permute.rare", nreps = 9)


#### Gamma richness trends ###############
# # permute rare
# gamma.trend.nat <- gamma.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.rare", nreps = 999)
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


#### Gamma above/below #######

# # Alpha differences above vs. below
#  system.time(alpha.above.trend <- div.part.alpha.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "permute", nreps =999))
# 
#  # Gamma differences with permute all : keeps beta diversities
#   system.time(gamma.above.nat.permute.all <- div.part.gamma.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "permute.all", nreps =999))
# 
# # Gamma differences with resample all : resamples beta diversity
#  system.time(gamma.above.nat.beta <- div.part.gamma.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "resample.beta", nreps =999))
# 
# save(alpha.above.trend,
#      gamma.above.nat.permute.all,gamma.above.nat.beta,
#      file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.covariable+year.Rdata")

# Used version:
 load (file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.withcofactors.Rdata") ##unimproved grasslands only

# Latest version: strangely, non significant corr btwn Nb.Acrit and Deltagamma.c
## Look into why it is different, probably a question of random null model => RERUN with more repetitions
# load(file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.covariable+year.Rdata") ## covariables including year of sampling

# Previous versions:
# load (file= "saved Rdata/article 2 - threshold/diversity.partitioning.Rdata") # all lucas grasslands
# load (file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.Rdata") # Unimproved without cofactors

# load (file= "saved Rdata/article 2 - threshold/diversity.partitioning.unimproved.with cofactors+bldg.Rdata") # unimproved grasslands + cofactors including building distance

#### gamma trends ALIEN richness ###############

#  # permute rare
#   gamma.trend.ali <- gamma.trend(spnames = rownames(glmSRali.overall$mean), group = aliens,
#                                  null.model="permute.rare", nreps = 999)
#   save(gamma.trend.ali,
#   file = "saved Rdata/article 2 - threshold/gamma.trends.SRali.Rdata")
#  load(file = "saved Rdata/article 2 - threshold/gamma.trends.SRali.Rdata")
#
#

#### Gamma and alpha above/below  ALIEN richness #######
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


#### Beta diversity = turnover + nestedness   => USED   ###########
 #  library(vegan)
 # betasor.nat = betasor.multi(spnames = impsp, group = natives, null.comm = T, bootstrap = F, nreps = 999)
 # betasor.boot.nat = betasor.multi(spnames = impsp, group = natives, null.comm = F, bootstrap = T, nreps = 999)
 
 
#save(betasor.nat,betasor.boot.nat, file = "saved Rdata/beta diversity output.3.Rdata") 
 load(file = "saved Rdata/beta diversity output.3.Rdata") #with covariables + year

# Old versions:
# load(file="saved Rdata/article 2 - threshold/betasor results.Rdata")
# load(file = "saved Rdata/beta diversity output.2.1.Rdata") ##with cofactors and RUN2

table.sorensen <-  as.data.frame(t(data.frame(lapply(betasor.nat, FUN = function(x) x[3,c("obs.below", "obs.above", "obs.diff", "z.diff", "p.diff")]))))
table.sorensen.boot <-  as.data.frame(t(data.frame(lapply(betasor.boot.nat, FUN = function(x) x[3,c("obs.below", "mean.below","sd.below", "p.below", "obs.above","mean.above","sd.above", "p.above","obs.diff", "mean.diff","sd.diff", "p.diff")]))))


table.turnover <- t(data.frame(lapply(betasor.nat, FUN = function(x) x[1,c("obs.below", "obs.above", "obs.diff", "z.diff", "p.diff")])))
table.nestedness <- as.data.frame(t(data.frame(lapply(betasor.nat, FUN = function(x) x[2,c("obs.below", "obs.above", "obs.diff", "z.diff", "p.diff")]))))
table.nestedness.boot <-  as.data.frame(t(data.frame(lapply(betasor.boot.nat, FUN = function(x) x[2,c("obs.below", "mean.below","sd.below", "p.below","obs.above","mean.above","sd.above", "p.above","obs.diff", "mean.diff","sd.diff", "p.diff")]))))


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




#### Extract result table for loss in native diversity above thresholds:   ######
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

# out$deltagamma.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$deltagamma[m,out[m, "th"]]  )
# out$delta.z.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$z.delta[m,out[m, "th"]]  )
# out$delta.null.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$null.delta[m,out[m, "th"]]  )
# out$GRP.beta <- sapply(sp.names, FUN = function(m)  gamma.above.nat.beta$P.delta[m,out[m, "th"]] )

out$BRo <- sapply(sp.names , FUN = function(m) gamma.above.nat.permute.all$betap.below[m,out[m, "th"]] )
out$BRc <-sapply(sp.names, FUN = function(m) gamma.above.nat.permute.all$betap.above[m,out[m, "th"]] )

out$z.beta.diff <-table.turnover[sp.names,"z.diff"]
out$P.beta.diff <-table.turnover[sp.names,"p.diff"]

# final formatting of table
table.div.part <- out[impsp,]


write.csv(as.matrix(orderBy(~ th.CI, table.div.part)), file="output diversity partitioning.csv", fileEncoding= "native.enc")


#### Extract result table for loss in alien diversity above thresholds:   #####
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







#### Table S1 with AR and NR critical abundances  ##############
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




#### TABLE 1 GLMS (only significant species) : 

# # OLD #####
# table1 <- data.frame( nobs = glmSRnat.overall$impact.spread[impsp,]$prevalence,
#                       glmSRnat.overall$glms[impsp,],
#                       DEM10.c = glmSRnat.overall$covar.tab$DEM_10[impsp,]$coef.glm,
#                       DEM10.P = glmSRnat.overall$covar.tab$DEM_10[impsp,]$P.coef,
#                       
#                       SLOPE.c = glmSRnat.overall$covar.tab$SLOPE[impsp,]$coef.glm,
#                       SLOPE.P = glmSRnat.overall$covar.tab$SLOPE[impsp,]$P.coef,
# 
#                       Northness.c = glmSRnat.overall$covar.tab$Northern[impsp,]$coef.glm,
#                       Northness.P = glmSRnat.overall$covar.tab$Northern[impsp,]$P.coef,
# 
#                       SRali.c = glmSRnat.overall$covar.tab$SRali[impsp,]$coef.glm,
#                       SRali.P = glmSRnat.overall$covar.tab$SRali[impsp,]$P.coef,
#                     
#                       abun.meanc = rowMeans(glmSRnat.overall$est[impsp,], na.rm = T)
# )

## formatted for JEcol in April 2018:
sel <- c(impsp, impsp.nocovar[!impsp.nocovar %in% impsp])

table1 <- data.frame( nocovar = glmSRnat.overall.nocovar$glms[sel,],
                      abun.meanc.nocovar = round(rowMeans(exp(glmSRnat.overall.nocovar$est[sel,])*100-100, na.rm = T), 1),
                      glmSRnat.overall$glms[sel,], 
                      Year.c = round(exp(glmSRnat.overall$covar.tab$year[sel,]$coef.glm)*100-100, 1),
                      Year.P = round(glmSRnat.overall$covar.tab$year[sel,]$P.coef, 4),
                      
                      DEM10.c = round(exp(glmSRnat.overall$covar.tab$DEM_10[sel,]$coef.glm)*100 -100 , 1),
                      DEM10.P = round(glmSRnat.overall$covar.tab$DEM_10[sel,]$P.coef,4),
                      
                      SLOPE.c = round(exp(glmSRnat.overall$covar.tab$SLOPE[sel,]$coef.glm)*100-100, 1),
                      SLOPE.P =round(glmSRnat.overall$covar.tab$SLOPE[sel,]$P.coef, 4),
                      
                      Northern.c = round(exp(glmSRnat.overall$covar.tab$Northern[sel,]$coef.glm)*100-100, 1),
                      Northern.P = round(glmSRnat.overall$covar.tab$Northern[sel,]$P.coef, 4),
                      
                      
                      SRali.c = round(exp(glmSRnat.overall$covar.tab$SRali[sel,]$coef.glm)*100-100, 1),
                      SRali.P = round(glmSRnat.overall$covar.tab$SRali[sel,]$P.coef, 4),
                      
                      abun.meanc = round(rowMeans(exp(glmSRnat.overall$est[sel,])*100-100, na.rm = T), 1)
)

table1$pdev.DEM_10 <- round((table1$dev.DEM_10 / table1$null.dev), 4)
table1$pdev.SLOPE <- round((table1$dev.SLOPE / table1$null.dev), 4)
table1$pdev.Northern <- round((table1$dev.Northern / table1$null.dev), 4)
table1$pdev.Year <- round((table1$dev.year / table1$null.dev), 4)
table1$pdev.SRali <- round((table1$dev.SRali / table1$null.dev), 4)
table1$pdev.abun <- round((table1$dev.abun / table1$null.dev), 4)


table1$dAICfull <- round(table1$aic.null - table1$aic.abun, 1)

table1$dAICabun <- round(table1$aic.SRali - table1$aic.abun, 1)

table1 <- table1[, c("nocovar.df", "nocovar.dev.ratio", "abun.meanc.nocovar",
        "df",  "dev.ratio", "dAICfull","dAICabun",
        "Year.c",  "pdev.Year", "DEM10.c",  "pdev.DEM_10", "SLOPE.c",  "pdev.SLOPE","Northern.c",  "pdev.Northern", "SRali.c",  "pdev.SRali","abun.meanc",  "pdev.abun")]

table1 <- table1[c(order(table.div.part$th.CI), 9,10,11),]
write.csv(table1, file= "table1.csv", row.names = row.names(table1))



# Table S1 : environmental drivers of SR

# Native Richness
f0<- glm( SRnat ~ 1, data = tmp, family = poisson)
f1<- glm( SRnat ~ DEM_10, data = tmp, family = poisson)
f2<- glm( SRnat ~ DEM_10 + SLOPE, data = tmp, family = poisson)
f3<- glm( SRnat ~ DEM_10 + SLOPE + Northern, data = tmp, family = poisson)
f4<- glm( SRnat ~ DEM_10 + SLOPE + Northern + SRali  , data = tmp, family = poisson)
f5<- glm( SRnat ~ DEM_10 + SLOPE + Northern +SRali + year, data = tmp, family = poisson)

AIC(f0,f1,f2,f3,f4, f5)
anova(f5)
summary(f5)

# Alien richness
f0<- glm( SRali ~ 1, data = tmp, family = poisson)
f1<- glm( SRali ~ DEM_10, data = tmp, family = poisson)
f2<- glm( SRali ~ DEM_10 + SLOPE, data = tmp, family = poisson)
f3<- glm( SRali ~ DEM_10 + SLOPE + Northern, data = tmp, family = poisson)
f4<- glm( SRali ~ DEM_10 + SLOPE + Northern + year, data = tmp, family = poisson)
f5<- glm( SRali ~ DEM_10 + SLOPE + Northern + year + SRnat, data = tmp, family = poisson)


AIC(f0,f1,f2,f3, f4, f5) 
anova(f4)
summary(f4) # best model has no interactions





#### TABLE S3 GLMS (ALL ALIEN focal species): ####
focal.aliens <- rownames(glmSRnat.overall$glms)
focal.aliens <- focal.aliens[focal.aliens %in% aliens]


tableS3init <- data.frame( nobs = glmSRnat.overall.nocovar$impact.spread[focal.aliens,]$prevalence ,
                       glmSRnat.overall.nocovar$glms[focal.aliens,],
                       abun.meanc = rowMeans(glmSRnat.overall.nocovar$est[focal.aliens,], na.rm = T),
                       abun.minP = apply(glmSRnat.overall.nocovar$P[focal.aliens,],1, min, na.rm = T)
)



tableS3 <- data.frame( nobs = glmSRnat.overall$impact.spread[focal.aliens,]$prevalence,
                      glmSRnat.overall$glms[focal.aliens,],
                      DEM10.c = glmSRnat.overall$covar.tab$DEM_10[focal.aliens,]$coef.glm,
                      DEM10.P = glmSRnat.overall$covar.tab$DEM_10[focal.aliens,]$P.coef,
                      
                      SLOPE.c = glmSRnat.overall$covar.tab$SLOPE[focal.aliens,]$coef.glm,
                      SLOPE.P = glmSRnat.overall$covar.tab$SLOPE[focal.aliens,]$P.coef,
                      
                      Northness.c = glmSRnat.overall$covar.tab$Northern[focal.aliens,]$coef.glm,
                      Northness.P = glmSRnat.overall$covar.tab$Northern[focal.aliens,]$P.coef,
                      
                      SRali.c = glmSRnat.overall$covar.tab$SRali[focal.aliens,]$coef.glm,
                      SRali.P = glmSRnat.overall$covar.tab$SRali[focal.aliens,]$P.coef,
                      
                      SRali.c = glmSRnat.overall$covar.tab$year[focal.aliens,]$coef.glm,
                      SRali.P = glmSRnat.overall$covar.tab$year[focal.aliens,]$P.coef,
                      
                      abun.meanc = exp(rowMeans(glmSRnat.overall$est[focal.aliens,], na.rm = T))
)

write.csv(cbind(tableS3init, tableS3), file= "tableS3.csv", row.names = row.names(tableS3))

### TABLES2 WITH LOGIT Function for focal species abundance correlation with covariables #####

library(MASS)
tableS2.focalsp <- data.frame(matrix(NA, length(focal.aliens), 16),row.names = focal.aliens)
colnames(tableS2.focalsp) <- c("species", "resid.df","%dev.best","DeltaAIC.best", "P.best", "best",
                               "elev.coef","slope.coef","north.coef","SRali.coef", "Year.coef",
                               "elev.P","slope.P","north.P","SRali.P","Year.P")

options(contrasts = c("contr.treatment", "contr.poly"))

for (i in 1:length(focal.aliens)) {
  sp <- focal.aliens[i] 
  abun <- tmp.comm[, sp]
  abun[abun == 0] <- NA
  elev <- tmp.envplot$DEM_10[!is.na(abun)]
  slope <- tmp.envplot$SLOPE[!is.na(abun)]
  SRali <- tmp.envplot$SRali[!is.na(abun)]
  SRnat <- tmp.envplot$SRnat[!is.na(abun)]
  north <- tmp.envplot$Northern[!is.na(abun)]
  Year <- tmp.envplot$year[!is.na(abun)]
  abun <- na.omit(abun)
  
  fit0 <- polr(as.factor(abun) ~ 1,  Hess=TRUE)
  fit5 <- polr(as.factor(abun) ~elev  + slope +  north + SRali + Year ,  Hess=TRUE)
  fit5.2 <- stepAIC(fit5, ~.^2, scope = list(upper = ~ elev  + slope +  north + SRali + Year, lower = ~1))
  AOV <- anova(fit0,fit5.2, fit5)
  n.var = length(fit5.2$coefficients)
  

  # Old version based on coeffs and Pvalues of the full model:
  # if (n.var == 5) AOV = anova(fit0, fit5)
  #coeffs
  # (ctable <- coef(summary(fit5)))
  #Odds ratio = exp(coef) for a logit function
  # OR <- exp(ctable[1:5, 1])
  # coef <- ctable[1:5, 1]
  ## calculate and store p values
  # p <- (pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2)[1:5]
  
  
  # Coeff and P values from the BEST model 5.2

  p <- as.vector(rep(1,5))
  names(p) = c("elev","slope","north","SRali","Year")
  if (n.var >0) { p [names(fit5.2$coefficients)] <- round(Anova(fit5.2)[,3],4) }
  
  coef <- as.vector(rep(0,5))
  names(coef) = c("elev","slope","north","SRali","Year")
  if (n.var >0) {  coef [names(fit5.2$coefficients)] <- round(fit5.2$coefficients,4) }
  
  # Illustrate effect of Year!
  par(mfrow = c(1,2))
  plot(sort(unique(Year)), tapply(abun, Year, mean))
  plot(slope, abun)
  
  #output
  tableS2.focalsp[i, ] <- c(species[sp, "SpeciesName"],
                            fit5.2$df.residual,(AOV[1,3] - AOV[2,3])/AOV[1,3] ,AIC(fit0) - AIC(fit5.2),
                            AOV[2,7],paste(unlist(fit5.2$terms), collapse = " "),
                            coef, p
                            )
  
  rm(fit0, fit5, fit5.2)
}


write.csv(tableS2.focalsp, file= "table S2 for focal species.csv")




#### TABLE 2 in the main text : ####
table2 <- data.frame(Acrit = table.div.part$th.CI,
                     Presence = table.div.part$prevalence,
                     above.Acrit = table.div.part$n.plot.impact,
                     dominance = table.div.part$n.plot.dominant,
                     alpha.below = round(table.div.part$aRo, 1),
                     alpha.below.sd = round(table.div.part$aRo.sd,1),
                     alpha.above = round(table.div.part$aRc,1),
                     alpha.above.sd = round(table.div.part$aRc.sd,1),
                     delta.alpha = table.div.part$aRc -  table.div.part$aRo,
                     perc.delta.alpha = round(((table.div.part$aRc -  table.div.part$aRo)/table.div.part$aRo)*100, 1),
                     gamma.below = table.div.part$GRo,
                     gamma.above = table.div.part$GRc,
                     delta.gamma.c = table.div.part$deltagamma - table.div.part$delta.null.permute.all,
                     delta.gamma.P = round(table.div.part$GRP.permute.all, 4),
                      delta.beta.z = table.turnover[rownames(table.div.part),"z.diff"],
                      delta.beta.P = table.turnover[rownames(table.div.part),"p.diff"]
                     
)
rownames(table2) <- table.div.part$species
table2 <- table2 [order(table.div.part$th.CI),]

write.csv(table2, file = "table2.csv")



# save results   ############
# save.image("saved Rdata/article 2 - threshold/article threshold 1.3.5.Rdata") # March 2018
# save.image("saved Rdata/article 2 - threshold/article threshold 1.3.6.Rdata") # April 9th 2018
# save.image("saved Rdata/article 2 - threshold/article threshold 1.3.7.Rdata") # July 2018
# save.image("saved Rdata/article 2 - threshold/article threshold 1.3.8.Rdata") # August 2018




### save bundle for publication:
nreps <- 999
db<- databp[databp$PlotName %in% unimprovedgrasslands,]
min.occur <- 5
min.class <- 2

glms_NR_output <- glmSRnat.overall

# data set to publish:
dataset = db[, c( "SurveyName",  "SurveyStartYear","SurveyEndYear",
                  "SpeciesCode","SpeciesName.x","NATIVE","ALIEN",
                  "PlotName","DominanceClassDescription", "abun",
                  "POINTX","POINTY","DEM_10","ASPECT","SLOPE", "Northern", "Eastern",
                  "SRnat", "SRali", "SR")]
a <- names(which( (rowSums(table(dataset$SpeciesCode, dataset$abun)[,2:6]>=min.occur)>=min.class)
                  &  table(dataset$SpeciesCode, dataset$abun)[,1]>=min.occur))
## Selecting sites where at least one of the focal species is present:
dataset <- dataset[which(dataset$SpeciesCode %in% a),]


save( dataset = dataset,
      glms_NR_output = glmSRnat.overall,
      boot.output = boot.output,
      boot.indices = boot.indices,
      file = "Bernard-Verdier&Hulme_JEcol2018_data.Rdata"
)

