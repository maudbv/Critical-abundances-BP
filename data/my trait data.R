# Import trait data from Banks Peninsula (from summer 2014/2015)
library(doBy)

# import Seed mass data (in g)
mySM <- read.csv(file="data/traits/SM measurements.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)
mySM.mean <- summaryBy(SM.mg. ~ Sp.code, data= mySM, na.rm= T)

# import height data (in cm)
myH<- read.csv(file="data/traits/Height measurements.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)
myH$Hmax = apply(cbind(myH$Hrep,myH$Hveg), 1, FUN = max, na.rm=T)
myH$Hmax [myH$Hmax== -Inf] = NA
myH.mean <- summaryBy(Hveg + Hrep + Hmax ~ spcode, data= myH, na.rm= T)


## leaf traits:
# import leaf area photo numbers
leaf.pictures<- read.csv(file="data/traits/leaf area pictures.csv", sep = c(","),
                          na.str=c("","NA"), as.is=T, stringsAsFactor=F)
leaf.pictures$nb.leaves.num <- as.numeric(leaf.pictures$nb.leaves.num)

#import leaf mass data
leaf.mass <- read.csv(file="data/traits/leaf mass measurements.csv",na.str=c("","NA"), sep = c(","), as.is=T, stringsAsFactor=F)


# import leaf area measurements from ImageJ
files <- paste("data/traits/measures",dir(path = "data/traits/measures"), sep ="/")
area.raw <- do.call("rbind",  lapply(files,
                                     function(x) read.delim(x, header = T,sep="\t",stringsAsFactors = FALSE)))
write.csv(area.raw,"area.raw.csv")

area.sum <- summaryBy(Area~ Label,area.raw, FUN = sum)

# match leaf pictures with leaf area measurements
areas <-sapply(leaf.pictures$photo. , FUN = function(x) {
  
  ind <- grep(x, area.sum$Label)
  area <- NA
  # if the number is also asn iphone picture
  if (length(ind) >1 & length(grep("IMG", area.sum$Label[ind])) ==1 )  ind <- ind[-grep("IMG", area.sum$Label[ind])]
  # if the number exists also as a picture from Mica's camera 
  if (length(ind) >1 & length(grep("DSC", area.sum$Label[ind])) ==1 )  ind <- ind[-grep("DSC", area.sum$Label[ind])]
  
  # if we have a normal unique index of the picture, extract the total area:
  if (length(ind) ==1)  area <- area.sum$Area.sum[ind]
  
  # if the number has some strange measured shapes added by ImageJ => sum the areas
  if (length(ind)>1 & length(grep(":", area.sum$Label[ind])) ==(length(ind) -1) ) {
    area <- sum(area.sum$Area.sum[ind])
    ind <-NA
  }
  
  # if the number has some strange measured shapes added by ImageJ => sum the areas
  if (length(ind)>1 & length(grep("_alternate", area.sum$Label[ind])) == 1)  {
    area <- sum(area.sum$Area.sum[ind])
    ind <-NA
  }
  
  # if the number has some strange measured shapes added by ImageJ => sum the areas
  if (length(ind)>1 & length(grep("-", area.sum$Label[ind])) == 1)  {
    area <- sum(area.sum$Area.sum[ind])
    ind <-NA
  }
  
  #if non of the above is true :
   if (length(ind)>1)  warning(paste(length(ind), "measures for", x))
  return(area)
})

leaf.pictures$area.total <- areas[match(leaf.pictures$photo., names(areas))]

write.csv(leaf.pictures, file = "leaf.pictures.csv")

# SUM area and mass per subrep (sum the subpictures)
subrep_area<- summaryBy(nb.leaves.num + area.total ~ spcode +obs + rep + subrep ,id= ~ name + type ,leaf.pictures, FUN = sum, na.rm =F)
subrep_mass<- summaryBy(as.numeric(dry.mass..mg.)~ name + obs + rep + subrep,  leaf.mass, FUN = sum )


# merge leaf mass with leaf pictures
mySLA <- merge(subrep_area,subrep_mass, by.x=c("obs","rep","subrep", "name"),by.y=c("obs","rep","subrep", "name"))
mySLA <- orderBy(~ name + obs + rep + subrep, mySLA) 

names(mySLA) <- c("obs","rep","subrep", "name","spcode", "nb.leaves","area.total","type", "dry.mass")                        

# merge large leaves that are spread over different subreps (ferns in particular)
# obs numbers : 
# 22, 23, 93          Pteridium esculentum
# 66                  Polystichum vestitum
# 10, 11, 9, 8, 7     Polystichum oculatum
obs_large_leaves <- c(22, 23,93, 66,10,11,9,8,7)
tmp <-  summaryBy(.~ name + obs, mySLA[mySLA$obs %in% obs_large_leaves,], id =~ spcode + rep + subrep + type, FUN = sum)
names(tmp) <- c("name","obs",  "nb.leaves","area.total", "dry.mass","spcode","rep","subrep","type")                        
mySLA <- mySLA[- which(mySLA$obs %in% obs_large_leaves),]
mySLA <- rbind(mySLA, tmp[,names(mySLA)])

### Remove erroneous data : 
  # Galium propinquum
  mySLA <- mySLA[-which(mySLA$spcode == "GALPRO"),]
  
  # Asplenium hookerianum : keep only rounded leaf samples
  mySLA  <- mySLA [-grep("elongated", mySLA $type),]
  
  # Carmichaelia australis : keep only photosynthetic stems, not the small leaves
  mySLA <-mySLA[-which(mySLA$spcode == "CARAUS" & mySLA$type == "leaves"),]


# Calculate SLA and LA
mySLA$leaf.area <- mySLA$area.total/mySLA$nb.leaves # leaf area in mm2
mySLA$leaf.dry.mass <- mySLA$dry.mass/mySLA$nb.leaves  # dry mass in mg
mySLA$sla <- mySLA$area.total/mySLA$dry.mass # sla in mm2/mg == m2/Kg

# Calculate mean leaf traits per species
mySLA.mean <- summaryBy(leaf.area + leaf.dry.mass + sla ~ name + spcode ,data= mySLA, FUN = c(mean, sd), na.rm=T)
mySLA.mean$nreps <- summaryBy(obs ~ name + spcode ,data= mySLA, FUN = length) [,3]
mySLA.mean$nleaves <- summaryBy(nb.leaves ~ name + spcode ,data= mySLA, FUN = sum, na.rm=T) [,3]

### merge mean traits into one dataframe 
mytraits <- merge(mySM.mean, myH.mean, by.x="Sp.code", by.y="spcode", all=T)
mytraits$logHmax <- log(mytraits$Hmax.mean)
mytraits$logSM <- log(mytraits$SM.mg..mean)

names(mytraits) = c("Sp.code" ,"SM","Hveg","Hrep","Hmax","logHmax","logSM")  

mytraits$SLA <- mySLA.mean[match(mytraits$Sp.code, mySLA.mean$spcode), "sla.mean"]
