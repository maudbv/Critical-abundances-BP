### Sub sampling the plots of grasslands for funcitonal trait analysis


### FILTER 1: landcover = grassland   + FILTER 2:  dominant = herbaceous   ##############

robgrasslands=as.character(envplot$PLOTID[grep("grassland", ignore.case=T, envplot$landcover)])

tmp=databp[databp$DominanceRank==1,]
grasslands=as.character(tmp[which(as.character(tmp$Growth.forms) %in% c('GR','HR')),"PlotName"])

realgrasslands <- robgrasslands [robgrasslands %in% grasslands]  ## 751 plots of grassland

#### grassplots:
tmp= colSums(comm[which(rownames(comm) %in% as.character(realgrasslands)),]>0, na.rm=T)
traitdata$grassland.occur=tmp[match(rownames(traitdata), names(tmp))]

sum(traitdata$grassland.occur>0, na.rm=T) ## 466 species

####  target species identification  ######## 
# freq >10 in the grasslands, abundance spans 4 classes

# Robin's species


#### (final) FILTER 3: target species is present



#### Consequent species and trait data list modification :




################### plot #sp vs. freq of grassland occur #################
y=sapply(seq(1:max(traitdata$grassland.occur, na.rm=T)), FUN=function(x) {
  y=sum(traitdata$grassland.occur==x, na.rm=T)
  return(y)
})
x=seq(1:max(traitdata$grassland.occur, na.rm=T))

par(mar=c(4,4,1,4), las=1)
plot(x, y, log="x", ann=F, type="h") 
par(new=T)
plot(x, ((N-cumsum(y))/N)*100,type="l", col="red", log="x", axes=F, ann=F)
axis(4)
abline(v=9, lty="dashed", col="grey")
abline(h=(N-cumsum(y)[9])/N*100, lty="dashed", col="grey")

mtext(1, text="Frequency (nb. plots)", line=2)
mtext(2, text="Nb. of species", line=2.5, las=0)
mtext(4, text="% total species", line=2.5, las=0)

######## what sp are we excluding with a threshold of at least 10 occurrences?###############

####  Lost species? ######
traitdata$lifeform=as.factor(traitdata$lifeform)
kept=traitdata[which(traitdata$grassland.occur>9) ,]
lost=traitdata[which(traitdata$grassland.occur<10 & traitdata$grassland.occur>0),]

table(traitdata$NATIVE[which(traitdata$grassland.occur>0)]) # 224/466 sp = 52 % natives in total sp pool
table(kept$NATIVE) # 40 % natives in total sp pool
table(lost$NATIVE )  # 58% natives in the lost pool

# no signif differences in lifeforms
chisq.test(table(as.factor(kept$lifeform)),table(as.factor(lost$lifeform ))) 
chisq.test(table(as.factor(traitdata$lifeform[which(traitdata$grassland.occur>0)])),table(as.factor(lost$lifeform ))) # no signif differences in lifeforms

par(mfrow=(c(1,2)), mar=c(8,4,2,2))
barplot(table(as.factor(traitdata$lifeform[which(traitdata$grassland.occur>0)])), las=2)  # 58% natives in the lost pool
barplot(table(as.factor(kept$lifeform )), las=2) 
barplot(table(as.factor(lost$lifeform )), las=2)  

### abundance classes

tmp=apply(comm[realgrasslands,lost$Sp.code], 2,FUN=function(x) x==6)
write.csv(databp[databp$PlotName %in% rownames(tmp)[which(rowSums(tmp)>0)],], file="lost dominants.csv")

which(colSums(tmp)>0)


###### which plots loose SR?

