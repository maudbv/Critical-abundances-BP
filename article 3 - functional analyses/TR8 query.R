###### Qurying the TR8 package
library(TR8)
# trait list
trlist <- c("canopy_height","leaf_dmc","seed_mass")
# species list
splist <- unique( unlist(sub("_", " ", species$tip)))
# extract trait data
tr8extract <- tr8(splist, trlist,gui_config = F)
tr8results <- as.data.frame(tr8extract@results)

sum(!is.na(tr8results$canopy_height))
sum(!is.na(tr8results$seed_mass))

#add to traitdata
rownames(tr8results)=sub(" ", "_",rownames(tr8results))
traitdata$H_ledatr8 <- tr8results[match( traitdata$tip,rownames(tr8results)), "canopy_height"]*100
traitdata$SM_ledatr8 <- tr8results[match( traitdata$tip,rownames(tr8results)), "seed_mass"]
traitdata$LDMC_ledatr8 <- tr8results[match( traitdata$tip,rownames(tr8results)), "leaf_dmc"]

# 
# plot(traitdata$H_ledatr8, traitdata$H_LEDA, xlim=c(0,150)) ## inconsistent Hcan !! ?
# abline(0,1)
# plot(traitdata$SM_ledatr8, traitdata$SM_LEDA)   ### consistent SM
# abline(0,1)
# 

