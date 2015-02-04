# Import Federico's Banks Peninsula Data
library(doBy)
importBPdata= function() {

species <- read.csv(file="data/Banks Peninsula_681_species_list.csv", sep=";",as.is=T, stringsAsFactors=F, fileEncoding="latin1")
ranks <- read.csv(file="data/Banks_Peninsula_community_data.csv", row.names="PLOTNAME",as.is=T, stringsAsFactors=F, fileEncoding="latin1" )
envplot <- read.csv(file="data/Banks_Peninsula_1227_plots_variables.csv",as.is=T, stringsAsFactors=F, fileEncoding="latin1")[,-1]

rownames(envplot) <- envplot$PLOTID
rownames(species) <- species$Sp.code
match(rownames(ranks), rownames(envplot))

######  creating a clear discimitating taxon nomenclature 
species$SpeciesName =species$Species
species$Genus=sapply(as.character(species$SpeciesName), FUN=function(x) {
  g=strsplit(x, split=" ")[[1]][1]
 return(g)})

species$Species=sapply(as.character(species$SpeciesName), FUN=function(x) {
  g=strsplit(x, split=" ")[[1]][2]
 return(g)})

# species names as in the tips of the phylgoeny
species$tip<-paste(species$Genus,species$Species, sep="_")
## differentiate tips for phylogeny (genus_species) from unique subspecies
species$tipcomplete <-apply(cbind(species$Genus,species$Species,species$subsp), 1 ,
                            function(x)  paste(x,sep="_", collapse="_"))

# species with different subspecies
duplic=rbind(species[which(duplicated(species[, c("Genus","Species")])),])
#               species[which(duplicated(species[, c("Genus","Species")], fromLast=TRUE)),])
duplic=orderBy(~ Genus + Species, duplic)

duplic$subspecies=sapply(as.character(duplic$SpeciesName), FUN=function(x) {
  g=strsplit(x, split=c(" ") )[[1]]
  g=paste(g[3:length(g)], collapse="-")
 return(g)})
species$subspecies='NA'
species[match(duplic$Sp.code,species$Sp.code),"subspecies"]=duplic$subspecies
species[match(duplic$Sp.code,species$Sp.code),"tip"]= paste(duplic$Genus,duplic$Species, sep="_")

#####  Subset of species : alien vs. native
tmp=species[species$ALIEN==1,]
aliens=as.character(tmp$Sp.code)

tmp=species[species$NATIVE==1,]
natives=as.character(tmp$Sp.code)

tmp=species[species$WEEDOC==1,]
weedocs=as.character(tmp$Sp.code)

#####  Species in the surveys for which no data is provided:
which(!names(ranks)%in%rownames(species))
names(ranks)[216]      # Only one species: ECHHIS = Echium hispanicum or/ hispidissimum or/ hispidum ?
ranks[,216]            # present in 1 plot at rank 33
ranks <- ranks[,-216]   # We can remove that species ?
ranks <- ranks[,rownames(species)]  # order the species the same way


##### transforming some functional groups
species$lifestyle=as.character(species$LIFE_DURATION)
species$lifestyle[species$lifestyle%in% c('A','A.B','A.P','B')]='A'
species$lifestyle[species$lifestyle%in% c('B.P','P')]='P'

species$woody=0
species$woody[species$Growth.forms %in% c("SH", "TR")]=1

########### Wilson's abundance surveys

survey <- read.csv(file="data/Banks Peninsula Botanical Survey 1983-88_Ranks data.csv")

# addind a column with the ordered dominance class levels
tmp=na.omit(unique(survey[,c('DominanceClass','DominanceClassDescription')]))
tmp=tmp[order(tmp[,1]),]
#                                              "domlevels"    "abun"
# 43               A                  Abundant 3              5
# 25             A-B           Abundant-Common 4              4
# 947            A-D         Abundant-Dominant 2              6
# 1                B                    Common 5              3
# 81             B-C         Common-Occasional 6              2
# 12               C                Occasional 7              1
# 675              D                  Dominant 1              6 
# 224            N/A        No Dominance Class 8              NA  
tmp$level=c(3,4,2,5,6,7,1,8)
survey$domlevels=survey$DominanceClass
levels(survey$domlevels) <- c(2,3,1,4,5,6,1,7)  # effectively merges levels 1 and 2
survey$domlevels=as.numeric(as.character(survey$domlevels))
survey$domlevels[survey$domlevels==7]=NaN


### Correcting species names

#   nomen <- read.csv(file="data/Nomenclatura_new.csv", as.is=T, stringsAsFactor=F)
#   names(nomen)=c("oldnames", "newnames")
#   nomen$old.SpeciesCode=sapply(strsplit(nomen$oldnames, split=" "), FUN=function(x) toupper(paste(substr(x[c(1)],1,3),substr(x[c(2)],1,3), sep=""))) 
#   nomen$new.SpeciesCode=sapply(strsplit(nomen$newnames, split=" "), FUN=function(x) toupper(paste(substr(x[c(1)],1,3),substr(x[c(2)],1,3), sep=""))) 
#   nomen$status = NA
#   nomen$status[nomen$old.SpeciesCode %in% survey$SpeciesCode]=1
# #   write.csv(nomen, file="nomen.csv")
#   
#   wilson.names=unique(survey[,c("SpeciesName", 'SpeciesCode', 'old.SpeciesCode')])
#   write.csv(wilson.names, file='wilson.names.csv')
  # replacing new species codes
#   survey$old.SpeciesCode = survey$SpeciesCode
#   survey$SpeciesCode = nomen[match(survey$SpeciesCode, nomen$old.SpeciesCode), "new.SpeciesCode"] 
#   survey$SpeciesCode[is.na(survey$SpeciesCode)] = as.character(survey$old.SpeciesCode[is.na(survey$SpeciesCode)])

  # UPDATED nomenclature, with additional changes by maud:
  nomen.add <- read.csv(file="data/nomen.additional.csv", as.is=T, stringsAsFactor=F)
  nomen.add[which(!nomen.add$wilson.code%in% survey$SpeciesCode),] # 5 unmatched species...
  nomen.add$status = NA
  nomen.add$status[nomen.add$wilson.code %in% survey$SpeciesCode]=1

  survey$SpeciesCode.tmp = survey$SpeciesCode
  survey$SpeciesCode = nomen.add[match(survey$SpeciesCode, nomen.add$wilson.code), "fede.code"] 
  survey$SpeciesCode[is.na(survey$SpeciesCode)] = as.character(survey$SpeciesCode.tmp[is.na(survey$SpeciesCode)])


# Comparing federico's ranks dataset with Hugh Wilson's
# 62 species have been added in federico's data (58 vascular, 4 non-vascular)
add.sp=species[which(!as.character(species$Sp.code)%in%as.character(survey$SpeciesCode)),]  # 62 species have been added in federico's data (58 vascular, 4 non-vascular)
write.csv(add.sp, file="add.sp.csv")

rm.sp=survey[which(!survey$SpeciesCode%in%names(ranks)),] 
rm.sp=rm.sp[-which(duplicated(rm.sp$SpeciesName)),] 
write.csv(rm.sp, file="rm.sp.csv")

table(rm.sp$SpeciesType) # 437 species have been removed; of which 291 non vasculars and 145 vascualrs
table(survey[-which(duplicated(survey$SpeciesName)),] $SpeciesType) # there was 295 non-vasulcars in the entire dataset !
nv.kept=survey[which(survey$SpeciesType=="Non-Vascular" & !survey$SpeciesName%in%rm.sp$SpeciesName),]


# Merge abundance with environmental data
m=merge(survey, envplot, by.x=c('PlotName'), by.y=c('PLOTID'))
databp=m
databp=merge(databp, species, by.x='SpeciesCode', by.y='Sp.code')

## REMOVE DUPLICATED lines corresponding to a combination of species x plot => where do these duplicate originate from ??
databp[which(duplicated(databp[, c('PlotName','SpeciesCode')])),]
databp=databp[which(!duplicated(databp[, c('PlotName','SpeciesCode')])),]

### UPDATING plot information 
envplot=orderBy(~PLOTID,envplot)
rownames(envplot)=envplot$PLOTID

# removing plots which have no data in survey (but they do in "ranks", thats strange no?)
old_envplot=envplot
envplot=envplot[sort(unique(as.character(databp$PlotName))),]

# Add species richness info to plot
envplot$SR=table(as.character(databp$PlotName))[as.character(envplot$PLOTID)]
envplot$SRali=tapply(databp$ALIEN,as.character(databp$PlotName), FUN=sum, na.rm=T)[as.character(envplot$PLOTID)]
envplot$SRnat=tapply(!databp$ALIEN,as.character(databp$PlotName), FUN=sum, na.rm=T)[as.character(envplot$PLOTID)]
envplot$SRwee=tapply(databp$WEEDOC,as.character(databp$PlotName), FUN=sum, na.rm=T)[as.character(envplot$PLOTID)]

envplot$SR.gr=tapply(databp$Growth.forms=='GR',as.character(databp$PlotName), FUN=sum, na.rm=T)[as.character(envplot$PLOTID)]
envplot$SR.hr=tapply(databp$Growth.forms=='HR',as.character(databp$PlotName), FUN=sum, na.rm=T)[as.character(envplot$PLOTID)]
envplot$SR.wood=tapply(databp$Growth.forms%in% c('TR','SH') ,as.character(databp$PlotName), FUN=sum, na.rm=T)[as.character(envplot$PLOTID)]
envplot$SR.fe=tapply(databp$Growth.forms=='FE',as.character(databp$PlotName), FUN=sum, na.rm=T)[as.character(envplot$PLOTID)]
envplot$SR.legumes=tapply(databp$Family=='Fabaceae',as.character(databp$PlotName), FUN=sum, na.rm=T)[as.character(envplot$PLOTID)]


## Defining Vegetation types
tmp=databp[databp$DominanceRank==1,]
grasslands=as.character(tmp[which(as.character(tmp$Growth.forms) %in% c('GR','HR')),"PlotName"])
woodlands=as.character(tmp[which(as.character(tmp$Growth.forms) %in% c('SH','TR')),"PlotName"])

lowlands=as.character(tmp[which(as.character(tmp$Growth.forms) %in% c('GR','HR') & tmp$DEM_10 < 400),"PlotName"])
highlands=as.character(tmp[which(as.character(tmp$Growth.forms) %in% c('GR','HR') & tmp$DEM_10 >= 400),"PlotName"])

## Alternative vegetation types:

# grasslands=as.character(envplot[which((envplot$SR -envplot$SR.wood)/envplot$SR >= 0.70),"PLOTID"])
# woodlands=as.character(envplot[which((envplot$SR -envplot$SR.wood)/envplot$SR<0.70),"PLOTID"])
# 
# 
# grasslands=as.character(envplot[which((envplot$SR -envplot$SR.wood)/envplot$SR >= 0.70),"PLOTID"])
# woodlands=as.character(envplot[which((envplot$SR -envplot$SR.wood)/envplot$SR<0.70),"PLOTID"])
# lowlands=as.character(envplot[which((envplot$SR -envplot$SR.wood)/envplot$SR >= 0.70 & envplot$DEM_10 < 400),"PLOTID"])
# highlands=as.character(envplot[which((envplot$SR -envplot$SR.wood)/envplot$SR < 0.70 & envplot$DEM_10 >= 400),"PLOTID"])


envplot$vegtype=NA
envplot$vegtype[as.character(envplot$PLOTID) %in% woodlands]='W'
envplot$vegtype[as.character(envplot$PLOTID) %in% grasslands]='G'

envplot$newX=ceiling(envplot$POINTX/1000)
envplot$newY=ceiling(envplot$POINTY/1000)


## import land-use classification by Robin :

# "La couche utilis?e est celle du Minist?re de l'Environnement bas?e sur une interpr?tation 
# d'images a?riennes datant de 1996-1997 (huit ans apr?s la fin de l'inventaire de H. Wilson)."

landcover <- read.csv(file="data/landcover.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)
envplot$landcover=landcover$LANDCOVER[match(envplot$PLOTID, landcover$PLOTID)]

## realgrasslands
tmp=databp[databp$DominanceRank==1,]
realgrasslands=as.character(tmp[which(as.character(tmp$Growth.forms) %in% c('GR','HR')) & which(as.character(tmp$landcover) %in% c('High Producing Exotic Grassland','Low Producing Grassland')),"PlotName"])


## REAL FINAL MERGIING
# Final merging
m=merge(survey, envplot, by.x=c('PlotName'), by.y=c('PLOTID'))
databp=m
databp=merge(databp, species, by.x='SpeciesCode', by.y='Sp.code')
databp$domlevels=as.numeric(as.character(databp$domlevels))
## REMOVE DUPLICATED lines corresponding to a combination of species x plot => where do these duplicate originate from ??
databp[which(duplicated(databp[, c('PlotName','SpeciesCode')])),]
databp=databp[which(!duplicated(databp[, c('PlotName','SpeciesCode')])),]

databp$abun=7-as.numeric(as.character(databp$domlevels))

## Community data matrix
tmp=databp
x=rep(unique(na.omit(as.character(tmp$PlotName))), 679)
x=sort(x)
y=rep(x=unique(as.character(tmp$SpeciesCode)),length(unique(na.omit(as.character(tmp$PlotName)))))

db=data.frame(PlotName=x, SpeciesCode=y)
s=merge(db, tmp[, c("PlotName", 'SpeciesCode', "abun")], by.x=c("PlotName", 'SpeciesCode'), by.y=c("PlotName", 'SpeciesCode'), all.x=T)

comm=t(unstack(s, form= abun ~ PlotName))
rownames(comm)=sort(unique(na.omit(as.character(tmp$PlotName))))
colnames(comm)=unique(as.character(tmp$SpeciesCode))
comm[is.na(comm) & !is.nan(comm)]=0  ### we keep the NAN information when species appear in the survey but have no abundance class (only 16 observations)
comm=as.data.frame(comm)

# presence-absence community data
# occur1=ranks
# occur1=ceiling(occur1/1000)

occur2=comm
occur2=ceiling(occur2/1000)


# list of objects to return
return(list(databp, comm, envplot, species, grasslands, woodlands,lowlands, highlands, occur=occur2, aliens=aliens, natives=natives, realgrasslands) )
}

d=importBPdata()
databp=d[[1]]
comm=d[[2]]
envplot=d[[3]]
species=d[[4]]
grasslands=d[[5]]
woodlands=d[[6]]
lowlands=d[[7]]
highlands=d[[8]]
occur=d[[9]] 
aliens=d[[10]]
natives=d[[11]]
realgrasslands=d[[12]]
rm(d, importBPdata)


### correct species names
# species$preferredtip <- names.species[match(species$tipcomplete,as.character(names.species$originalname)),"preferredname"]
# species$preferredtip <- sapply(species$preferredname,FUN=function(x) paste(strsplit(x, " " )[[1]][c(1,2)],collapse="_" ))
# species$preferredtip[species$preferredtip %in% c("NA_NA", "")] <- species$tip[species$preferredtip %in% c("NA_NA", "")]

