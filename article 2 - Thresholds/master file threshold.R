## Threshold master file

#setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R")

library(doBy)
library(vegan)


# load("saved Rdata/article 2 - threshold/article threshold 1.0.Rdata")

### import data
# load("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/saved Rdata/article 3 - trait and phylo/save article 3.Rdata.RData")

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


### threshold analysis using GLMs

# set bootstrap sample size
R <- 99
db <- databp[databp$PlotName %in% realgrasslands,]
min.occur <- 5
min.class <- 2

#####  Bootstrpping within each class :

#  system.time(glmSR <- glm.test(db = db,bootstrap = T, min.occur= min.occur, R=R, CI=0.95, drastic =F))
#  glmSRnat <- glm.test(db = db,var="SRnat",bootstrap = T,min.occur= min.occur,  R=R, CI=0.95, drastic =F)
#  glmSRali <- glm.test(db = db,var="SRali",bootstrap = T,min.occur= min.occur,  R=R, CI=0.95, drastic =F)

# save(glmSR,glmSRnat,glmSRali, file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")
load(file = "saved Rdata/article 2 - threshold/booststrapped.glms.Rdata")

#####  Overall bootstrapping :

# create bootstrapped samples
  boot.output <- bootstrap.dataset(db=db, min.occur =min.occur,  min.class = min.class, R = R)
  save( boot.output , file = "saved Rdata/article 2 - threshold/boot.output.2classes.Rdata")

# load(file = "saved Rdata/article 2 - threshold/boot.output.Rdata")

## calculate glms on bootstraps
system.time(glmSR.overall <- glm.overallboot(db = db,boot.output=boot.output, variable = 'SR', min.occur= min.occur, min.class = min.class, R=R))

system.time(glmSRnat.overall <- glm.overallboot(db = db,boot.output=boot.output, variable = 'SRnat', min.occur= min.occur, min.class = min.class, R=R))
system.time(glmSRali.overall <- glm.overallboot(db = db,boot.output=boot.output, variable = 'SRali', min.occur= min.occur, min.class = min.class, R=R))

save(glmSR.overall, glmSRnat.overall, glmSRali.overall, file = "saved Rdata/article 2 - threshold/overall.boot.glms.2classe.Rdata")

# load(file = "saved Rdata/article 2 - threshold/overall.boot.glms2.Rdata")


# ## GLMs with dominance index as covariate
# system.time(glmSR.dom <- glm.test(db = db,bootstrap = T, R=R,min.occur= min.occur,  covar ="dominance"))
# glmSRnat.dom <- glm.test(db = db,var="SRnat",bootstrap = T, R=R, min.occur= min.occur, covar ="dominance")
# glmSRali.dom <- glm.test(db = db,var="SRali",bootstrap = T, R=R,min.occur= min.occur,  covar ="dominance")
# save(glmSR.dom,glmSRnat.dom,glmSRali.dom, file = "saved Rdata/article 2 - threshold/booststrapped+dominance.glms.Rdata")


## overall bootstrap impact size 

impact.SR <- impact.size (glmSR.overall)

impact.SRnat <- impact.size (glmSR.overall)

impact.SRali <- impact.size (glmSR.overall)

# ### Transform results by restricting target species :
# db <- databp[databp$PlotName %in% realgrasslands,]
# min.occur <- 5
# min.class <- 2
# a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=min.class) 
#                   &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
# 
# glmSR.overall   <- lapply (glmSR.overall , FUN = function(X) X[a,])
# 
# glmSRnat.overall  <- lapply (glmSRnat.overall, FUN = function(X) X[a,])
# 
# glmSRali.overall <- lapply (glmSRali.overall, FUN = function(X) X[a,])


### SUMMARY of frequencies of significantly negative effects and thresholds
glmSR.sum <-  summary.glmtest(M=glmSR.overall, group="ALIEN", type="overall.boot")
glmSRnat.sum <- summary.glmtest(M=glmSRnat.overall, group="ALIEN", type="overall.boot") 
glmSRali.sum <- summary.glmtest(M=glmSRali.overall, group="ALIEN", type="overall.boot") 




######## calculate proportional impacts (mean % of species lost)
glmSR$boot.thresh = add.prop(N = glmSR, var ="SR", data=db)
glmSRnat$boot.thresh = add.prop(N = glmSRnat, var ="SRnat", data=db)
glmSRali$boot.thresh = add.prop(N = glmSRali, var ="SRali", data=db)

# save results
 save.image("saved Rdata/article 2 - threshold/article threshold 1.1.Rdata")

