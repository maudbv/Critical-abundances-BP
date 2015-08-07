## Threshold master file

setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R_alienimpactBP/")
# setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R")


# load packages
library(doBy)
library(vegan)
library(markdown)

### import data
# load("saved Rdata/article 2 - threshold/article threshold 1.2.Rdata")

source('script/data/import BP species and environment data.R', encoding = "native.enc")
source('script/data/import trait data.R', encoding = "native.enc")

### taxonomy solving :
# extract names, match them to DB, create a reference list of changes, update initial DB


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



### threshold analysis using GLMs

# set bootstrap sample size
nreps <- 999
db <- databp[databp$PlotName %in% realgrasslands,]
min.occur <- 5
min.class <- 2

#####  Bootstrpping within each class :

#  system.time(glmSR <- glm.test(db = db,bootstrap = T, min.occur= min.occur, nreps=nreps, CI=0.95, drastic =F,min.occur =min.occur,  min.class = min.occur))
#  glmSRnat <- glm.test(db = db,var="SRnat",bootstrap = T,min.occur= min.occur,  nreps=nreps, CI=0.95, drastic =F, min.occur =min.occur,  min.class = min.occur)
#  glmSRali <- glm.test(db = db,var="SRali",bootstrap = T,min.occur= min.occur,  nreps=nreps, CI=0.95, drastic =F,min.occur =min.occur,  min.class = min.occur)

# save(glmSR,glmSRnat,glmSRali, file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")
#load(file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")

#####  Overall bootstrapping :

## create bootstrapped samples
# system.time(boot.output <- bootstrap.dataset(db=db, min.occur =min.occur,  min.class = min.class, nreps = nreps))
# system.time(boot.indices <- extract.indices(boot.output, db = db))
#
# save( boot.output, boot.indices ,
#  file = "saved Rdata/article 2 - threshold/boot.output.2.0.Rdata")

load(file = "saved Rdata/article 2 - threshold/boot.output.2.0.Rdata")

# save( boot.output, boot.indices ,
#         file = "saved Rdata/article 2 - threshold/boot.output.99reps.Rdata")

# calculate glms on bootstraps
# system.time(glmSR.overall <- glm.overallboot(db = db,boot.indices=boot.indices, sp.target = sp.target, variable = 'SR', min.occur= min.occur, min.class = min.class, nreps=nreps))
# # system.time(glmSRnat.overall <- glm.overallboot(db = db,boot.indices=boot.indices, sp.target = sp.target, variable = 'SRnat', min.occur= min.occur, min.class = min.class, nreps=nreps))
#  system.time(glmSRali.overall <- glm.overallboot(db = db,boot.indices=boot.indices, sp.target = sp.target, variable = 'SRali', min.occur= min.occur, min.class = min.class, nreps=nreps))
#
# save(glmSR.overall, glmSRnat.overall, glmSRali.overall,file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.0.Rdata")

  load(file = "saved Rdata/article 2 - threshold/overall.boot.glms.2.0.Rdata")


# ## GLMs with dominance index as covariate
# system.time(glmSR.dom <- glm.test(db = db,bootstrap = T, nreps=nreps,min.occur= min.occur,  covar ="dominance"))
# glmSRnat.dom <- glm.test(db = db,var="SRnat",bootstrap = T, nreps=nreps, min.occur= min.occur, covar ="dominance")
# glmSRali.dom <- glm.test(db = db,var="SRali",bootstrap = T, nreps=nreps,min.occur= min.occur,  covar ="dominance")
# save(glmSR.dom,glmSRnat.dom,glmSRali.dom, file = "saved Rdata/article 2 - threshold/booststrapped+dominance.glms.Rdata")


# ### Transform results by restricting target species :
# db <- databp[databp$PlotName %in% realgrasslands,]
# min.occur <- 5
# min.class <- 2
# a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=min.class)
#                   &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
#
# glmSR.overall   <- lapply (glmSR.overall , FUN = function(X) X[a,])
#
# glmSRali.overall <- lapply (glmSRali.overall, FUN = function(X) X[a,])
#
#
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

######## calculate proportional impacts (mean % of species lost)
glmSR$boot.thresh = add.prop(N = glmSR, var ="SR", data=db)
glmSRnat$boot.thresh = add.prop(N = glmSRnat, var ="SRnat", data=db)
glmSRali$boot.thresh = add.prop(N = glmSRali, var ="SRali", data=db)



 # Select pecies which show a threshold for Native richness
 impsp <- rownames(glmSRnat.overall$impact.spread[which(!is.na(glmSRnat.overall$impact.spread$th.CI)
                                                        & (rownames(glmSRnat.overall$impact.spread) %in% aliens)),])

## gamma and beta trends
source('script/article 2 - Thresholds/gamma.trend.R')
alpha.trend.nat <- alpha.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.rare", nreps = 999)
gamma.trend.nat <- gamma.trend(spnames = rownames(glmSRnat.overall$mean), null.model="permute.rare", nreps = 999)
beta.trend.nat <- beta.trend(spnames = rownames(glmSRnat.overall$mean),null.model="permute.rare", nreps = 999)
save(gamma.trend.nat, beta.trend.nat, file = "saved Rdata/article 2 - threshold/gamma.trends.Rdata")

load(file = "saved Rdata/article 2 - threshold/gamma.trends.Rdata")

# ##### Partition Gamma/beta/alpha native diversity for each imp species above/below #######
source('script/article 2 - Thresholds/div partitioning.R')
system.time(gamma.above.trend <- div.part.gamma.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "permute", nreps =999))
system.time(alpha.above.trend <- div.part.alpha.nm(spnames = rownames(glmSRnat.overall$mean), group=natives, null.model = "permute", nreps =999))

#system.time(divpart.ali.perm <- div.part.nm(spnames = impsp, group=aliens, null.model = "permute", nreps = 999))

save(gamma.above.trend, alpha.above.trend,file= "saved Rdata/article 2 - threshold/diversity.partitioning.Rdata")

load (file= "saved Rdata/article 2 - threshold/diversity.partitioning.Rdata")

#### Beta diversity = turnover + nestedness
betasor.nat = betasor(spnames = impsp, group = natives, nreps = 499)
betasor.ali = betasor(spnames = impsp, group = aliens, nreps = 499)


#### SPATIAL CLUSTERING
# source('C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R - alienimpactBP/script/article 2 - Thresholds/spatial correlation.R')
# save(spDst, spatial.mnnd, mnnd.th, mnnd.trends, mnnd.above,
#      file = "saved Rdata/article 2 - threshold/spatial MNND.Rdata"))

load(file = "saved Rdata/article 2 - threshold/spatial MNND.Rdata")



# extract result table for diversity partitionning :
out <- glmSRnat.overall$impact.spread
out <- out[impsp,c("th.CI", "prevalence", "n.plot.impact")]
out<- cbind(species = species[impsp, "tip"], out)

out$aRo <- sapply(impsp, FUN = function(m) alpha.above.trend$alpha.below[m,out[m, "th.CI"]] )
out$aRc <- sapply(impsp, FUN = function(m) alpha.above.trend$alpha.above[m,out[m, "th.CI"]] )
out$aRnull <- sapply(impsp, FUN = function(m) alpha.above.trend$null.above[m,out[m, "th.CI"]]  )
out$aRsdnull <- sapply(impsp, FUN = function(m)  alpha.above.trend$sdnull.above[m,out[m, "th.CI"]]  )
out$alpha.loss <- out$aRnull  - out$aRc
out$aR.P <- sapply(impsp, FUN = function(m)  alpha.above.trend$P.above[m,out[m, "th.CI"]] )

out$GRo <- sapply(impsp, FUN = function(m) gamma.above.trend$gamma.below[m,out[m, "th.CI"]] )
out$GRc <- sapply(impsp, FUN = function(m) gamma.above.trend$gamma.above[m,out[m, "th.CI"]] )
out$GRnull <- sapply(impsp, FUN = function(m)  gamma.above.trend$null.above[m,out[m, "th.CI"]]  )
out$GRsdnull <- sapply(impsp, FUN = function(m)  gamma.above.trend$sdnull.above[m,out[m, "th.CI"]]  )
out$gamma.loss <- out$GRnull  - out$GRc
out$GR.P <- sapply(impsp, FUN = function(m)  gamma.above.trend$P.above[m,out[m, "th.CI"]] )

out$BRo <- sapply(impsp, FUN = function(m) gamma.above.trend$betap.below[m,out[m, "th.CI"]] )
out$BRc <-sapply(impsp, FUN = function(m) gamma.above.trend$betap.above[m,out[m, "th.CI"]] )
table.div.part <- out

table.div.part <- orderBy(~ th.CI, table.div.part)

write.csv(data.frame(table.div.part), file="output diversity partitioning.csv")

# save results
 save.image("saved Rdata/article 2 - threshold/article threshold 1.2.Rdata")

