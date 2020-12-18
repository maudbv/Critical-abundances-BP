#### LEDA trait data importation and correction

require(doBy)

# SLA data
ledaSLA     <- read.csv(file="data/LEDA SLA.csv", stringsAsFactor=F)
ledaSLA.mean=summaryBy(mean.SLA..mm.2.mg. ~ SBS.name, data=ledaSLA, FUN=c(mean, sd), na.rm=T)

# SM data
ledaSM      <- read.csv(file="data/LEDA SM.csv", stringsAsFactor=F)
ledaSM.mean <- summaryBy(mean.SM..mg. ~ SBS.name, data=ledaSM, FUN=c(mean, sd), na.rm=T)


# Canopy Height data
ledaHcan      <- read.csv(file="data/LEDA H.canopy.csv", stringsAsFactor=F)
ledaHcan.mean <- summaryBy(mean.CH..m. ~ SBS.name, data=ledaHcan, FUN=c(mean, sd), na.rm=T)

# merge three trait info
LEDAtraits <- merge(merge(ledaSLA.mean, ledaSM.mean, by="SBS.name"), ledaHcan.mean ,by="SBS.name")
LEDAtraits$tip <- sapply(as.character(LEDAtraits$SBS.name), FUN=function(x) {paste(strsplit(x, split=".")[[1]][1:2], collapse=".")})
LEDAtraits <- summaryBy(.~ tip, data=LEDAtraits,id="SBS.name",FUN=c(mean), na.rm=T)
names(LEDAtraits) <-  c("SBS.name", "SLA.mean", "SLA.sd", "SM.mean", "SM.sd", "Hcan.mean", "Hcan.sd","tip")

LEDAtraits$match.tip <- 0 
LEDAtraits$match.tip[which(LEDAtraits$tip %in% species$tip)] <- 1

#### correct names to match our taxonomy
# names.species=resolve.names(df=LEDAtraits, colname="tip")

# 
# grep("Hordeum",LEDAtraits$tip)
# ledaSLA$tip[grep("Gnaphalium",ledaSLA$tip)]


###
write.csv(LEDAtraits, file="LEDA traits.csv")


