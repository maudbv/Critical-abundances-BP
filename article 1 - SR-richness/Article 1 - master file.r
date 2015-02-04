### Article 1 : abundance-richness relationships

setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R")

library(doBy)
library(vegan)

# import data
source('script/data/import fede data.R')

source('script/data/import trait data.R')

# save.image(file="saved Rdata/Banks peninsula dataset.Rdata")

# load functions
source('script/functions/p2star.R')
source('script/functions/function SR-class.R')
source('script/functions/simulation1.2.R')
source('script/functions/simulation1.1.R')
source('script/functions/simulation 2.R')
# source('script/functions/simulation3.R')
source('script/functions/combination.test.R')

# calculate raw relationships
source('script/article 1 - SR-richness/calculate class-SR relationships.R')

save(effect.summary, grass.effect.summary, wood.effect.summary,file="saved Rdata/article 1/effect summary.Rdata")

# null model testing
# source('script/article 1 - SR-richness/null model testing.R')
sapply(dir(path="saved Rdata/article 1/null model results"), FUN=function(x) load(file=paste("saved Rdata/article 1/null model results", x, sep="/"), envir = .GlobalEnv))

# SAVE data
# save.image(file="saved Rdata/JEcol resubmission.Rdata")
# save.image(file="saved Rdata/JEcol alternative vegtypes.Rdata")

# Plot graphs for JEcol article : see other script