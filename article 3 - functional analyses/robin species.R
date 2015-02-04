##### Robin's species list ########

robinsp <- read.csv(file="data/robinsp.csv",na.str=c("","NA"), as.is=T, sep=";" ,stringsAsFactor=F)


###### MATCHING with my species   ###############
robinsp$Sp.code <- unlist(sapply(robinsp$Species, FUN=function(x) {
 n= agrep(x, species$SpeciesName,max.distance=list(all=1), value=F)
if (!length(n)==1) {
   y <- NA
 }
else y <- as.character(species$Sp.code[n])
return(y)
}))

# 3 species without an obvious match :
dim(robinsp)
tmp <- robinsp$Species[which(!robinsp$Sp.code %in% species$Sp.code)]

# explore with fuzzy matching :
lapply(tmp, function(x) species[agrep(x, species$SpeciesName, max.distance=list(all=8), ignore.case=T, value=F), c("SpeciesName", "Sp.code")])

# "critesion murinum" = OK
robinsp$Sp.code[robinsp$Species == "Critesion murinum"] <- "CRIMUR"

# need more research for the two others :
# "Gnaphalium audax"   = Euchiton audax (preferred synonym)
## http://www.nzflora.info/factsheet/Taxon/Euchiton_audax.html
robinsp$Sp.code[robinsp$Species == "Gnaphalium audax"] <- as.character(species$Sp.code[species$tip=="Euchiton_audax"])

### "Leontodon taraxacoides" = synonym of Leontodon saxatilis
# http://www.ars-grin.gov/cgi-bin/npgs/html/taxon.pl?454196
robinsp$Sp.code[robinsp$Species == "Leontodon taraxacoides"] <- as.character(species$Sp.code[species$tip=="Leontodon_saxatilis"])

robinsp$tip <-  species$tip[match(robinsp$Sp.code, species$Sp.code)]
row.names(robinsp) <- robinsp$Sp.code

###### COUNTING how many species we have traits for ########

dim(robinsp[which( robinsp$Sp.code%in%rownames(traitdata)[traitdata$SLA == 1] ),]) # 65
dim(robinsp[which( robinsp$Sp.code%in%rownames(traitdata)[traitdata$Hmax == 1] ),]) # 81
dim(robinsp[which( robinsp$Sp.code%in%rownames(traitdata)[traitdata$SM == 1] ),]) # 63


### information on robins species occurrence
tmp= unique(as.character(databp[which((databp$tip %in% robinsp$tip) & (databp$PlotName %in% grasslands) ),"PlotName"]))
tmp=colSums(comm[tmp,])
traitdata$robin.grassoccur=tmp[match(rownames(traitdata), names(tmp))] 
## = 537 species in both grasslands and robin sp occupied sites


# 21 robin species not on the phylogeny
robinsp[which(! robinsp$tip %in% BPtree$tip),"Species"]





