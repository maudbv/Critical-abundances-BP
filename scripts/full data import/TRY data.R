### TRY database query and formatting

# # preparing query sheet for trait data
TRYsp<-read.csv(file="data/Try species list.csv", stringsAsFactor=F)
TRYsp$Genus     <- sapply(TRYsp$AccSpeciesName, function(x) strsplit(x, split=" ")[[1]][1])

# 
# TRYsp$tip <- sapply(as.character(TRYsp$AccSpeciesName), FUN=function(x) {
# paste(strsplit(x, split=" ")[[1]][1:2], collapse="_")})
# length(which( species$tip %in% unique(TRYsp$tip)))                  #  478 of species
# table(species[which( species$tip %in% unique(TRYsp$tip)),"ALIEN"])  #  215 natives
# 
# TRYsp$match=0
# TRYsp$match[which(TRYsp$tip %in% species$tip)]=1
# write(paste(TRYsp[TRYsp$match==1,"AccSpeciesID"][1:200], collapse=","), file="1_200accessIDlist.txt")
# write(paste(TRYsp[TRYsp$match==1,"AccSpeciesID"][201:400], collapse=","), file="201_400accessIDlist.txt")
# write(paste(TRYsp[TRYsp$match==1,"AccSpeciesID"][401:481], collapse=","), file="401_481accessIDlist.txt")


# REsults from query of all trait data
TRY1 <- read.csv(file="data/Try1.csv", stringsAsFactor=F, row.names="Species")
TRY2 <- read.csv(file="data/Try2.csv", stringsAsFactor=F, row.names="Species")
TRY3 <- read.csv(file="data/Try3.csv", stringsAsFactor=F, row.names="Species")

trait.names <- sort(unique(c(names(TRY1), names(TRY2), names(TRY3))))
sp.tip      <- TRYsp[TRYsp$match==1,"tip"]


TRY <- data.frame(matrix(NA, nrow=length(sp.tip), ncol=length(trait.names)), row.names=sp.tip)
colnames(TRY)<- trait.names

TRY[rownames(TRY1), names(TRY1)]<- TRY1
TRY[rownames(TRY2), names(TRY2)]<- TRY2
TRY[rownames(TRY3), names(TRY3)]<- TRY3
TRY[is.na(TRY)]<- 0

# Potential species lists:
TRYsp.SLA   <- rownames(TRY)[which(TRY$Specific_leaf_area_.SLA.>0)]
TRYsp.Hveg  <- rownames(TRY)[which(TRY$Plant_height_vegetative>0)]
TRYsp.SM   <- rownames(TRY)[which(TRY$Seed_mass>0)]


# # EXPLORE:
# 
# length(which(TRY$Seed_mass>0))  # 362 SM
# length(which(TRY$Specific_leaf_area_.SLA.>0)) #333
# length(which(TRY$Plant_height_vegetative>0)) #375
# length(which(TRY$Plant_height_generative>0)) #196
# 
# 
# # Potential LHS trait dataset with TRY :
# TRYsp.match   <- rownames(TRY)[which(TRY$Plant_height_vegetative*TRY$Specific_leaf_area_.SLA.*TRY$Seed_mass>0)]
# #260 species with all three traits
# 
# TRYsp.match   <- rownames(TRY)[which(TRY$Plant_height_vegetative + TRY$Specific_leaf_area_.SLA.+ TRY$Seed_mass>0)]
# #441 species with at least one trait
# 
# table(species[species$tip %in% TRYsp.match, "ALIEN"])  # 186 Natives

# TRYsp.match[!(TRYsp.match %in% traitdata$tip[!is.na(traitdata$SLA)])]
# sum(!is.na(traitdata$SLA))
# TRYsp.match[!(TRYsp.match %in% traitdata$tip[!is.na(traitdata$Hmax)])]
# sum(!is.na(traitdata$Hmax))
# TRYsp.match[!(TRYsp.match %in% traitdata$tip[!is.na(traitdata$SM)])]
# sum(!is.na(traitdata$SLA))

# ## Potential genus list added to dataset:
# 
# TRYgen.sla=unique(TRYsp$Genus[which(! (TRYsp$Genus[match(TRYsp.SLA, TRYsp$tip)]) %in% rownames(gentrait)[gentrait$SLA>0])])
# 
# sum(!(unique(TRYsp$Genus[match(TRYsp.SLA, TRYsp$tip)]) %in% rownames(gentrait)[gentrait$SLA>0]))
# #85
# 
# sum(!(unique(TRYsp$Genus[match(TRYsp.Hveg, TRYsp$tip)]) %in% rownames(gentrait)[gentrait$Hmax>0]))
# #71
# sum(!(unique(TRYsp$Genus[match(TRYsp.SM, TRYsp$tip)]) %in% rownames(gentrait)[gentrait$SM>0]))
# #64