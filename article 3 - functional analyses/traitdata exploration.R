### Exploring trait data for Banks peninsula
load("saved Rdata/article 1/effect summary.Rdata")


db=databp[databp$PlotName %in% realgrasslands,]
a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=min.class) 
                  &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))


N=dim(traitdata[which(traitdata$grassland.occur>15 & rownames(traitdata)%in% a ),])[1]
# 466

dim(traitdata[which(traitdata$grassland.occur>9),])[1]
# 153

chisq.test(table(as.factor(lost$lifeform )),table(as.factor(kept$lifeform )))


### change in sp richness 

envplot$SRsubset[match(row.names(comm), envplot$PLOTID)]=rowSums(comm[,which(!names(comm) %in% lost$Sp.code)]>0)

plot(SRsubset ~ SR, envplot)
abline(0,1)

# nb lost species per plot per landcover 
par(mar=c(4,16,2,1))
boxplot( envplot$SR- envplot$SRsubset ~ envplot$landcover, las=2, horizontal=T, xlab="Lost species richness")

# % lost species richness
par(mar=c(4,16,2,1))
boxplot( (envplot$SR- envplot$SRsubset)/(envplot$SR) * 100~ envplot$landcover, las=1, horizontal=T, xlab="% Lost species richness")


# % lost species richness
par(mar=c(4,16,2,1))
boxplot( (envplot$SR- envplot$SRsubset)/(envplot$SR) * 100~ envplot$vegtype, las=1, horizontal=T, xlab="% Lost species richness")




###  rarefaction for traits!   ####################
gaptolerance=function(th, TRY=F){
  n=NULL
  m=NULL
  s=NULL
  h=NULL
  sm=NULL
  p=NULL
  for (i in 0:th){
    
    if (TRY==F) {
    n=c(n,dim(traitdata[which(traitdata$grassland.occur>i),])[1])
    miss=traitdata[which(traitdata$grassland.occur>i & (traitdata$SLA + traitdata$Hmax + traitdata$SM)<3),]
    m=c(m,dim(miss)[1])
    h=c(h,sum(miss$Hmax==0))
    s=c(s,sum(miss$SLA==0))
    sm=c(sm,sum(miss$SM==0))
    p= (n-m)/n *100
    }
    
    if (TRY==T) {
      n=c(n,dim(traitdata[which(traitdata$grassland.occur>i),])[1])
      miss=traitdata[which(traitdata$grassland.occur>i & (traitdata$SLA.try + traitdata$Hmax.try + traitdata$SM.try)<3),]
      m=c(m,dim(miss)[1])
      h=c(h,sum(miss$Hmax.try==0))
      s=c(s,sum(miss$SLA.try==0))
      sm=c(sm,sum(miss$SM.try==0))
      p= (n-m)/n *100
    }
    
  }

  par(mar=c(4,4,1,4))
  plot(1,1,xlim=c(0,th), ylim=c(1,100), type="n", ann=F)
  title(ylab="% remaining species", xlab = "minimum frequency threshold")
  # lines(0:th, n, lty="dotted")
  lines(0:th,(n)/N*100, lty="dotted")
  
  par(new=T)
  plot(0:th, p, axes=F, ann=F, type="l", ylim=c(0, 100))
  axis(4)
  mtext(4,text="% species covered", las=0, line=2)
  
  lines(0:th, (n-s)/n*100, col="forestgreen")
  lines(0:th, (n-h)/n*100, col="firebrick")
  lines(0:th, (n-sm)/n*100, col="slateblue")
  
  legend("bottomleft", legend=c("remaining species", "total missing species", "missing SLA", "missing Hmax", "missing SM"),
         col=c("black", "black", "forestgreen","firebrick","slateblue"),
         lty=c("dashed","solid","solid","solid","solid"), cex=0.7, bty="n")
  
  grid()
}

gaptolerance(th=40)
gaptolerance(th=40, TRY=T)

##### SP OK species with all 3 traits:  ############################

spok=traitdata[which(traitdata$grassland.occur>0 & (traitdata$SLA + traitdata$Hmax + traitdata$SM)==3),]
### 205 species with all 3 trait values

sum(spok$ALIEN==0)  # Only  49 native species have all 3...
table(spok$Growth.forms)  # mostly herbs and some grasses
sum(rownames(spok) %in% rownames(effect.summary))  # 132 of the targets
dim(spok)[1]  

# mean and variation per genus
# gentrait= as.data.frame(sapply(c("SLA", "Hmax", "SM"),FUN=function(x) {
#   y=tapply(traitdata[,x],as.factor(traitdata$Genus), FUN=mean, na.rm=T)
#   return(y)
# }))


##### Tables of missing trait information for grassland species :

missing=traitdata[which(traitdata$grassland.occur>0 & (traitdata$SLA + traitdata$Hmax + traitdata$SM)<3),]
dim(missing)[1]
### 277 missing species in grasslands

table(missing$ALIEN) # 261 natives issing
sum(missing$Hmax==0) # 208heights missin
sum(missing$SLA==0) # 312 SLA missing
sum(missing$SM==0) # 286 SM missing

# compared to TRY potential ?
sum(TRYsp.SLA %in% missing$tip[which(missing$SLA==0)]) # 101/355
sum(TRYsp.Hveg %in% missing$tip[which(missing$Hmax==0)]) # 85/240 
sum(TRYsp.SM %in% missing$tip[which(missing$SM==0)]) # 92/302

misstry=traitdata[which(traitdata$grassland.occur>0 
                        & (traitdata$SLA.try + traitdata$Hmax.try + traitdata$SM.try)<3),]
# 288 entirely missing species
table(misstry$ALIEN) # 232 natives issing
sum(misstry$Hmax.try==0) # 151 heights missing
sum(misstry$SLA.try==0) # 248 SLA missing
sum(misstry$SM.try==0) # 219 SLA missing



### Only FREQUENT species :
dim(traitdata[which(traitdata$grassland.occur>15),])
# 153 sp

sum(rownames(spok) %in% rownames(traitdata[which(traitdata$grassland.occur>9),]))
# 109 LHS complete from the frequents
sum(rownames(spok.try) %in% rownames(traitdata[which(traitdata$grassland.occur>9),]))
# 103 LHS complete with TRY from the frequents

missing.imp=traitdata[which(traitdata$grassland.occur>9 & (traitdata$SLA + traitdata$Hmax + traitdata$SM)<3),]
### 44 frequent missing species 
table(missing.imp$ALIEN) 
sum(missing.imp$Hmax==0) 
sum(missing.imp$SLA==0) 
sum(missing.imp$SM==0) 

### TRY additions potentially
write.csv(missing.imp, "grassland species missing traits.csv")

