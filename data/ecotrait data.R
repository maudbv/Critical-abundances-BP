### importing ECOTRAIT database
### trait dataset curated by Margaret Watts for Landcare Research
### Sp.code should match those in Hugh Wilson's abundance surveys.

ecotraits <- read.csv(file="data/Ecotrait database.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)


## sp.code matching and verification:
if (length(which(!species$Sp.code %in% ecotraits$Sp.code))==0) {
  ## all the species in the survey are included in ecotraits
message("All species codes are included in the ecotrait dataset")
}
