### calculations for revisions JEcol

### proportion of weeds with inmpact
tmp=effect.summary[effect.summary$ALIEN==1,]
tmp=grass.effect.summary[grass.effect.summary$ALIEN==1,]
tmp$WEEDOC=species[rownames(tmp),"WEEDOC"]
table(tmp$WEEDOC,tmp$effect)
summary(table(tmp$WEEDOC,tmp$effect))
# Chisq = 0.6927, df = 2, p-value = 0.7073  => no effect

table(tmp$WEEDOC,tmp$nat.effect)
summary(table(tmp$WEEDOC,tmp$nat.effect))
# Chisq = 0.4615, df = 2, p-value = 0.7939

table(tmp$WEEDOC,ceiling(tmp$total.effect>0))
summary(table(tmp$WEEDOC,ceiling(tmp$total.effect>0)))
# Chisq = 2.4038, df = 2, p-value = 0.3006

table(tmp$WEEDOC,tmp$ali.effect)
summary(table(tmp$WEEDOC,tmp$ali.effect))




#### Vegetation types description




### woodlands vs grasslands
library(doBy)
dat=databp
dat5 = dat[which(dat$DominanceRank<5),]
dat1 = dat[which(dat$DominanceRank==1),]
ftable( dat1$vegtype, dat1$abun)
dim(dat1[which(dat1$abun<4 & dat1$vegtype=="W"),]) ## 91 plot where woody 1st ranked is not more than common
weakW=dat1[which(dat1$abun<4 & dat1$vegtype=="W"),"PlotName"]

# abundance class of the first ranked species:
x=t(as.matrix(table(dat1$abun, dat1$vegtype)))/colSums(table(dat1$abun, dat1$vegtype))
cumsum(x[1,6:1])
cumsum(x[2,6:1])



v=cbind(as.data.frame(ftable(dat$PlotName, dat$woody) ), as.data.frame(ftable(dat5$PlotName, dat5$woody) )[,3])
# v=as.data.frame(ftable(dat$PlotName[dat$DominanceRank<3], dat$woody[dat$DominanceRank<3] ) )
names(v)=c("PlotName", "woody","freq","freq5")
v$domclas=dat1[match(v$PlotName, dat1$PlotName),"abun"]
v$sp=dat1[match(v$PlotName, dat1$PlotName),"SpeciesName.y"]
v$SR=envplot$SR[match(v$PlotName, envplot$PLOTID)]
v$vegtype=envplot$vegtype[match(v$PlotName, envplot$PLOTID)]
v=orderBy(~ vegtype + PlotName + woody, v)
v$perc=v$freq/v$SR *100

boxplot(perc~  woody + vegtype, data=v)
(boxplot(freq~  woody + vegtype, data=v))
(boxplot(freq5~  woody + vegtype, data=v))
ftable( v$vegtype,v$woody, v$freq5)

# among woodlands, which plots have over 30% of woody species?
vw=v[which(v$vegtype=="W" & v$woody==1 ),]
vg=v[which(v$vegtype=="G" & v$woody==1 ),]
hist(vw$perc)
hist(vg$perc)

hist(vw$freq5)
hist(vg$freq5)

dim(vw)
tmp=vw[vw$perc>20,]
dim(tmp) # 260 ont plus de 10% d'especes woody 
tmp=vw[vw$perc<20,]
dim(tmp) # 37 grasslands ont plus de 10% d'especes woody 

# 91 woodlands have a less than abundant-common first ranking woody sp ( 80% have a first rank sp which is at least abundant-common)
# 91 + 59 have either an abundant first sp, or more than 30% of woody sp
weirdwood=vw[which(vw$PlotName %in% weakW & vw$perc<20 ),] # 59 of these have less than 30% of woody sp
weirdwood=vw[which(vw$PlotName %in% weakW & vw$perc<20 & vw$freq5<3 ),] # 48 have less than3 of the first ranking sp that are woody
dim(weirdwood)[1]/300 # 10% sont bizarrement classes
write.csv(weirdwood, file="weird woodlands.csv")
