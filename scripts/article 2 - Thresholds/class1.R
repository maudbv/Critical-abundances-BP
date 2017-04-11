######### Testing what happens in class 1
library(doBy)
library(lme4)

# Clear decrease in richness when aliens have equivalent classes of abundaces :
## HYP : low abundances of alien are usually also associated to high abundance/or high richness or other aliens ?

db=databp[which(!is.na(databp$abun) & !is.na(databp$ALIEN) & databp$abun==1),]
boxplot( SR ~  ALIEN , db)
boxplot( SR ~  vegtype , db)
boxplot( SR ~  ALIEN  + woody, db)
boxplot( SRnat~  ALIEN + woody + vegtype, db) ### only native trees reduce SR !!
boxplot( SRali ~  ALIEN + woody, db) ### only native trees reduce SRali

# Only aliens
db=databp[which(!is.na(databp$abun) &!is.na(databp$vegtype) & databp$ALIEN==1 & databp$abun==1),]

boxplot( SR ~  vegtype , db)
boxplot( SR ~  woody, db)
boxplot( SR ~  Growth.forms, db)
boxplot( SR ~   woody + vegtype, db)
boxplot( SRnat ~   woody + vegtype, db)
boxplot( SRali ~  woody + vegtype, db)

boxplot( SR ~  woody + vegtype, db)
summary(fm0 <- lmer( SR ~  (1|SpeciesCode), db))
summary(fm1 <- lmer( SR ~  vegtype + (1|SpeciesCode), db)) 
summary(fm2 <- lmer( SR~  woody + vegtype + (1|SpeciesCode), db)) 
AIC(fm0, fm1, fm2) ## seems a bit more significant...

boxplot( SRnat ~  woody + vegtype, db)
summary(fm0 <- lmer( SRnat ~  (1|SpeciesCode), db))
summary(fm1 <- lmer( SRnat ~  vegtype + (1|SpeciesCode), db)) 
summary(fm2 <- lmer( SRnat ~  woody + vegtype + (1|SpeciesCode), db)) 
AIC(fm0, fm1,  fm2) ## woodiness does not explain more

boxplot( SRali ~  woody + vegtype, db)
summary(fm0 <- lmer( SRali ~  (1|SpeciesCode), db))
summary(fm1 <- lmer( SRali ~  vegtype + (1|SpeciesCode), db)) 
summary(fm2 <- lmer( SRali ~  woody + vegtype + (1|SpeciesCode), db)) 
AIC(fm0, fm1,  fm2) ## woodiness does not explain much more


# mean values per species
db=databp[which(!is.na(databp$abun) & databp$ALIEN==1 & databp$abun==1),]
db.mean=summaryBy(SR + SRnat +SRali ~ SpeciesCode, id=~Growth.forms, data=db)
boxplot( SR.mean ~  Growth.forms, db.mean)
boxplot( SRnat.mean ~  Growth.forms, db.mean)
boxplot( SRali.mean ~  Growth.forms, db.mean)

# mixed model ?
summary(fm0 <- lmer( SR ~   (1|SpeciesCode), db))
summary(fm1a <- lmer( SR ~  as.numeric(abun) + (1|SpeciesCode) , db))
summary(fm1b <- lmer( SR ~  ALIEN + (1|SpeciesCode) , db))
summary(fm2 <- lmer( SR ~  as.numeric(abun) + ALIEN + (1|SpeciesCode), db))
summary(fm3 <- lmer( SR ~  as.numeric(abun) + ALIEN + Growth.forms + (1|SpeciesCode), db))
summary(fm4 <- lmer( SR ~  as.numeric(abun) + ALIEN + vegtype + (1|SpeciesCode), db))

AIC(fm0, fm1a, fm1b ,fm2, fm3, fm4)



boxplot( SRali ~  ALIEN + abun, databp)
summary(lm( SR ~  ALIEN*abun, databp))

AIC(

tmp=databp[databp$abun ==1,]

boxplot( SR ~  Growth.forms + ALIEN, tmp)
boxplot( SR ~  ALIEN, tmp)

class1.mean = data.frame(SR.mean=rep(NA, length(unique( tmp$SpeciesCode))))
class1.mean$SR.mean = tapply(tmp[, c("SR")], tmp$SpeciesCode, FUN="mean", na.rm=T)