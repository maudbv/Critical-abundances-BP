###### Importing TRait datasets
f=function(X=species, TRY=FALSE) {
require(doBy)

traitdata <- X

##############  NZ Ecotraits database
source('script/data/ecotrait data.R')
traitdata$PlantHeight_ecotrait  <- ecotraits[match(traitdata$Sp.code, ecotraits$Sp.code), "MeanPlantHeight"] * 100   # convert to cm
traitdata$Naturalisation        <- ecotraits[match(traitdata$Sp.code, ecotraits$Sp.code), "NaturalisationPeriod"]
traitdata$Naturalisation [traitdata$ALIEN == 0] <- NA
traitdata$Nat_start  <- as.numeric(as.vector(lapply(strsplit(traitdata$Naturalisation, split=" - "), FUN= function(x) x[1])))
traitdata$Nat_end    <- as.numeric(lapply(strsplit(traitdata$Naturalisation, split=" - "), FUN= function(x) x[2]))

traitdata$SeedLength_ecotrait <- ecotraits[match(traitdata$Sp.code, ecotraits$Sp.code), "SeedLength"]
traitdata$SM_ecotrait <- ecotraits[match(traitdata$Sp.code, ecotraits$Sp.code), "MeanSeedMass"]
traitdata$SM_ecotrait[traitdata$SM_ecotrait == 0] <- NA
traitdata$LeafLength_ecotrait <- ecotraits[match(traitdata$Sp.code, ecotraits$Sp.code), "MeanMaximumLeafLength"]
traitdata$LeafWidth_ecotrait <- ecotraits[match(traitdata$Sp.code, ecotraits$Sp.code), "LeafWidth"]
traitdata$LeafSize_ecotrait <- traitdata$LeafLength * traitdata$LeafWidth


lifeforms <- read.csv(file="data/lifeform correction list.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)
traitdata$lifeform <- lifeforms[match(traitdata$Sp.code, lifeforms$Sp.code), "lifeform_corrected"]

# write.csv(traitdata, file="corrected Ecotrait database")

# ####### Import LEDA trait data
# 
# SLA data
ledaSLA     <- read.csv(file="data/LEDA SLA.csv", stringsAsFactor=F)
ledaSLA$tip <- sapply(as.character(ledaSLA$SBS_name), FUN=function(x) {paste(strsplit(x, split="_")[[1]][1:2], collapse="_")})
ledaSLA_mean=summaryBy(mean_SLA_.mm.2.mg. ~ tip, data=ledaSLA, FUN=c(mean, sd), na.rm=T)
 traitdata$SLA_LEDA <- ledaSLA_mean[match( traitdata$tip,ledaSLA_mean$tip), "mean_SLA_.mm.2.mg..mean"]
              
# SM data
ledaSM      <- read.csv(file="data/LEDA SM.csv", stringsAsFactor=F)
ledaSM$tip  <- sapply(as.character(ledaSM$SBS_name), FUN=function(x) {paste(strsplit(x, split="_")[[1]][1:2], collapse="_")})
ledaSM_mean <- summaryBy(mean_SM_.mg. ~ tip, data=ledaSM, FUN=c(mean, sd), na.rm=T)
 traitdata$SM_LEDA <- round(ledaSM_mean[match( traitdata$tip,ledaSM_mean$tip), "mean_SM_.mg..mean"],4)
                    
# Canopy Height data
ledaHcan      <- read.csv(file="data/LEDA Hcanopy.csv", stringsAsFactor=F)
ledaHcan$tip  <- sapply(as.character(ledaHcan$SBS_name), FUN=function(x) { paste(strsplit(x, split="_")[[1]][1:2], collapse="_")})
ledaHcan_mean <- summaryBy(mean_CH_.m. ~ tip, data=ledaHcan, FUN=c(mean, sd), na.rm=T)
 traitdata$H_LEDA <- ledaHcan_mean[match( traitdata$tip,ledaHcan_mean$tip), "mean_CH_.m..mean"]*100
  

####### Nico Gross trait data (only Height and SLA)
NGross.data       <- read.csv(file="data/NGross trait data.csv", stringsAsFactor=F)
NGross.data$tip<- paste(NGross.data$Genus, NGross.data$species, sep="_")

NGross.data_mean <- summaryBy(VegetativeHeight + SLA ~ Species.Code, id=~tip,data=NGross.data , FUN=c(mean, sd), na.rm=T)
traitdata$H_NG <- NGross.data_mean[match( traitdata$Sp.code,NGross.data_mean$Species.code), "VegetativeHeight.mean"]
traitdata$H_NG[is.na(traitdata$H_NG)] <- NGross.data_mean[match( traitdata$tip,NGross.data_mean$tip), "VegetativeHeight.mean"][is.na(traitdata$H_NG)]
traitdata$H_NG[is.na(traitdata$H_NG)] <- NGross.data_mean[match( traitdata$preferredtip,NGross.data_mean$tip), "VegetativeHeight.mean"][is.na(traitdata$H_NG)]


traitdata$SLA_NG <- NGross.data_mean[match(  traitdata$Sp.code,NGross.data_mean$Species.code), "SLA.mean"]/10
traitdata$SLA_NG[is.na(traitdata$SLA_NG)] <- NGross.data_mean[match( traitdata$tip,NGross.data_mean$tip), "SLA.mean"][is.na(traitdata$SLA_NG)] /10
traitdata$SLA_NG[is.na(traitdata$SLA_NG)] <- NGross.data_mean[match( traitdata$preferredtip,NGross.data_mean$tip), "SLA.mean"][is.na(traitdata$SLA_NG)] /10



###  Ordonez 2014 trait data
source('script/data/ordonez trait data.R')
traitdata$SLA_ordo=ordonez[match(traitdata$preferredtip, ordonez$tip), "SLA"]
traitdata$SLA_ordo[is.na(traitdata$SLA_ordo)]=ordonez[match(traitdata$tip, ordonez$tip), "SLA"][is.na(traitdata$SLA_ordo)]

traitdata$Hmax_ordo=ordonez[match(traitdata$preferredtip, ordonez$tip), "Hmax"]
traitdata$Hmax_ordo[is.na(traitdata$Hmax_ordo)]=ordonez[match(traitdata$tip, ordonez$tip), "Hmax"][is.na(traitdata$Hmax_ordo)]

traitdata$SM_ordo=ordonez[match(traitdata$preferredtip, ordonez$tip), "SM"]
traitdata$SM_ordo[is.na(traitdata$SM_ordo)]=ordonez[match(traitdata$tip, ordonez$tip), "SM"][is.na(traitdata$SM_ordo)]



# #### Merging continuous data
# 
# traitdata$SLA=apply( cbind(traitdata$SLA_NG, traitdata$SLA_LEDA),1, FUN=mean, na.rm=T)
# traitdata$SLA.sd=apply( cbind(traitdata$SLA_NG, traitdata$SLA_LEDA),1, FUN=sd, na.rm=T)
# 
# #200 sp
# 
# traitdata$Hmax=apply( cbind(traitdata$H_NG, traitdata$H_LEDA,traitdata$PlantHeight_ecotrait),1, FUN=mean, na.rm=T)
# traitdata$Hmax.sd=apply( cbind(traitdata$H_NG, traitdata$H_LEDA,traitdata$PlantHeight_ecotrait),1, FUN=sd, na.rm=T)
# #344 sp
# 
# traitdata$SM=apply( cbind(traitdata$SM_LEDA,traitdata$SM_ecotrait),1, FUN=mean, na.rm=T)
# traitdata$SM.sd=apply( cbind(traitdata$SM_LEDA,traitdata$SM_ecotrait),1, FUN=sd, na.rm=T)
# #269 sp



### source TRY trait info
source('script/data/TRY data.R')

# ADD INFORMATION of TRY IN TRAIT DATA
traitdata$TRYsla<-NA
traitdata$TRYsla[match(TRYsp.SLA, traitdata$tip)] <- 1

traitdata$TRYhveg<-NA
traitdata$TRYhveg[match(TRYsp.Hveg, traitdata$tip)] <- 1

traitdata$TRYsm<-NA
traitdata$TRYsm[match(TRYsp.SM, traitdata$tip)] <- 1




# summarize available data:
traitdata$SLA= apply(!is.na(cbind(traitdata$SLA_NG,
                                  traitdata$SLA_LEDA,
                                  traitdata$SLA_ordo)),1,
                     function(x) as.numeric(sum(x)>0))
traitdata$Hmax= apply(!is.na(cbind(traitdata$H_NG,
                                   traitdata$H_LEDA,
                                   traitdata$Hmax_ordo,
                                   traitdata$PlantHeight_ecotrait)),1,
                      function(x) as.numeric(sum(x)>0))
traitdata$SM= apply(!is.na(cbind( traitdata$SM_LEDA,
                                  traitdata$SM_ordo,
                                  traitdata$SM_ecotrait)),1,
                    function(x) as.numeric(sum(x)>0))

traitdata$SLA.try= apply(cbind(traitdata$SLA,traitdata$TRYsla),1,
                         function(x) as.numeric(sum(x, na.rm=T)>0) )
traitdata$Hmax.try= apply(cbind(traitdata$Hmax,traitdata$TRYhveg),1,
                          function(x) as.numeric(sum(x, na.rm=T)>0) )
traitdata$SM.try= apply(cbind(traitdata$SM,traitdata$TRYsm),1,
                        function(x) as.numeric(sum(x, na.rm=T)>0) )





return(traitdata)
}

traitdata=f()
rm(f)

#additional modifications of traitdata
traitdata$ALIEN=species[rownames(traitdata),"ALIEN"]


# write.csv(traitdata, file="traitdata.csv")
