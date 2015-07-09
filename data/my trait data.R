# Import trait data from Banks Peninsula (from summer 2014/2015)

setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
#setwd("~/Dropbox/Work/doc boulot/post doc Lincoln/R")

# import Seed mass data
mySM <- read.csv(file="data/traits/SM measurements.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)
mySM.mean <- summaryBy(SM.mg. ~ Sp.code, data= mySM, na.rm= T)

# import height data
myH<- read.csv(file="data/traits/Height measurements.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)
myH$Hmax = apply(cbind(myH$Hrep,myH$Hveg), 1, FUN = max, na.rm=T)
myH$Hmax [myH$Hmax== -Inf] = NA
myH.mean <- summaryBy(Hveg + Hrep + Hmax ~ spcode, data= myH, na.rm= T)


# merge into one dataframe
mytraits <- merge(mySM.mean, myH.mean, by.x="Sp.code", by.y="spcode", all=T)
mytraits$logHmax <- log(mytraits$Hmax.mean)
mytraits$logSM <- log(mytraits$SM.mg..mean)

names(mytraits) = c("Sp.code" ,"SM","Hveg","Hrep","Hmax","logHmax","logSM")  

# compare to database values:
plot(logHmax~ logSM,data= mytraits, bg = "black", pch = 21)
abline(lm(logHmax~ logSM,data= mytraits),  bg = "black", pch = 21)
points(log(traitdata$H_LEDA)~ log(traitdata$SM_LEDA))
abline(lm(log(traitdata$H_LEDA)~ log(traitdata$SM_LEDA)), lty = "dotted")

# import leaf area photo numbers
leaf.pictures<- read.csv(file="data/traits/leaf area pictures.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)
leaf.mass <- read.csv(file="data/traits/leaf mass measurements.csv",na.str=c("","NA"), as.is=T, stringsAsFactor=F)

tmp <- merge(leaf.pictures, leaf.mass, by.x=c("obs","rep", "name"),by.y=c("obs","rep", "name"))
tmp <- orderBy(~ obs + rep + subrep + subpicture, tmp) 



tmp <- merge(myH, leaf.pictures, by.x=c("obs","rep", "spcode", "name", "date"),by.y=c("obs","rep","spcode", "name","date"))
tmp <- orderBy(~ obs + rep + subrep + subpicture, tmp) 
      