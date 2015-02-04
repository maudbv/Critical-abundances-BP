#### LEDA trait data importation and correction

require(doBy)

# SLA data
ledaSLA     <- read.csv(file="data/LEDA SLA.csv", stringsAsFactor=F)
ledaSLA_mean=summaryBy(mean_SLA_.mm.2.mg. ~ SBS_name, data=ledaSLA, FUN=c(mean, sd), na.rm=T)

# SM data
ledaSM      <- read.csv(file="data/LEDA SM.csv", stringsAsFactor=F)
ledaSM_mean <- summaryBy(mean_SM_.mg. ~ SBS_name, data=ledaSM, FUN=c(mean, sd), na.rm=T)


# Canopy Height data
ledaHcan      <- read.csv(file="data/LEDA Hcanopy.csv", stringsAsFactor=F)
ledaHcan_mean <- summaryBy(mean_CH_.m. ~ SBS_name, data=ledaHcan, FUN=c(mean, sd), na.rm=T)

# merge three trait info
LEDAtraits <- merge(merge(ledaSLA_mean, ledaSM_mean, by="SBS_name"), ledaHcan_mean ,by="SBS_name")
LEDAtraits$tip <- sapply(as.character(LEDAtraits$SBS_name), FUN=function(x) {paste(strsplit(x, split="_")[[1]][1:2], collapse="_")})
LEDAtraits <- summaryBy(.~ tip, data=LEDAtraits,id="SBS_name",FUN=c(mean), na.rm=T)
names(LEDAtraits) <-  c("SBS_name", "SLA.mean", "SLA.sd", "SM.mean", "SM.sd", "Hcan.mean", "Hcan.sd","tip")

LEDAtraits$match_tip <- 0 
LEDAtraits$match_tip[which(LEDAtraits$tip %in% species$tip)] <- 1

#### correct names to match our taxonomy
# names.species=resolve.names(df=LEDAtraits, colname="tip")

# 
# grep("Hordeum",LEDAtraits$tip)
# ledaSLA$tip[grep("Gnaphalium",ledaSLA$tip)]


###
write.csv(LEDAtraits, file="LEDA traits.csv")


