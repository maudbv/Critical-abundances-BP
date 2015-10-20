### import trait data from Ordonez et al. 2014

ordonez <- read.csv(file="data/traits/ordonez2014 supplementary file.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)
dim(ordonez)

# # pattern matching with our species names
# tmp1=lapply(ordonez$Species,FUN=function(i) agrep(i, species$SpeciesName,
#                                                  max.distance = list(all=0)))
# sum(unlist(lapply(tmp1, length))) #296
# 
# match.ordonez=data.frame(rbind(cbind(which(unlist(lapply(tmp1, length))==1), 
#       unlist(tmp1[which(unlist(lapply(tmp1, length))==1)]),
#       match=1),
# cbind(which(unlist(lapply(tmp1, length))==2), 
#       unlist(lapply(tmp1[which(unlist(lapply(tmp1, length))==2)], paste, collapse="_")),
#       match=2)))
# names(match.ordonez)=c("ordonez.names", "mynames")
# 
# # with 1 letter fuzzy
# tmp2=lapply(ordonez$Species,FUN=function(i) agrep(i, species$SpeciesName,
#                                                  max.distance = list(all=1)))
# sum(unlist(lapply(tmp2, length))) #308
# 
# 
# fmatch.ordonez=data.frame(rbind(cbind(which(unlist(lapply(tmp2, length))==1), 
#                                      unlist(tmp2[which(unlist(lapply(tmp2, length))==1)]),
#                                      match=1),
#                                cbind(which(unlist(lapply(tmp2, length))==2), 
#                                      unlist(lapply(tmp2[which(unlist(lapply(tmp2, length))==2)], paste, collapse="_")),
#                                      match=2)))
# names(fmatch.ordonez)=c("ordonez.names", "mynames")
# 
# # fuzzy errors:
# fuz=ordonez$Species[which(!fmatch.ordonez$ordonez.names %in% match.ordonez$ordonez.names)]
# 
# lapply(fuz,FUN=function(i) agrep(i, species$SpeciesName, value=T,max.distance = list(all=3)))
# 
# ###### Does not work well... and not very interesting in any case, only 12 additional species with fuzzy matching


### Add tips
ordonez$tip <- sub(" ", "_",ordonez$Species)

# trait averaging and simplification
ordonez$SLA=apply(ordonez[,c("SLA.Native","SLA.Alien")], 1, mean, na.rm=T)
ordonez$Hmax=apply(ordonez[,c("Hmax.Native","Hmax.Alien")], 1, mean, na.rm=T)
ordonez$SM=apply(ordonez[,c("SWT.Native","SWT.Alien")], 1, mean, na.rm=T)

# 
# 
# ### comparisons
# ordonez[which(ordonez$tip%in% missing$tip),]  # 83 species match our missing species
# dim(ordonez[which(ordonez$tip%in% species$tip),]) # 291 matching "tip", but 296 using "grep"
# 
# length(na.omit(ordonez[which(ordonez$tip%in% missing$tip),"SLA"])) # 71 missing SLA values
# length(na.omit(ordonez[which(ordonez$tip%in% missing$tip),"Hmax"])) # 82 missing H values
# length(na.omit(ordonez[which(ordonez$tip%in% missing$tip),"SM"])) # 65 missing SM values
