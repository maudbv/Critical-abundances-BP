# Species richness - species dominance class relationships

## Testing the relationships for each species with kruskal wallis + spearman
classSR=classSR.fun(db=databp, var='SR',domclass='domlevels',zeros=T,envplot=envplot, min.occur=10, min.class=3, alpha=0.01)
classSRnat=classSR.fun(db=databp, zeros=T, var='SRnat',min.occur=10, envplot=envplot,min.class=3, alpha=0.01)
classSRali=classSR.fun(db=databp, zeros=T, var='SRali', min.occur=10, envplot=envplot,min.class=3, alpha=0.01)
classSRwee=classSR.fun(db=databp, zeros=T, var='SRwee', min.occur=10,envplot=envplot, min.class=3, alpha=0.01)

classSR.grass=classSR.fun(db=databp[databp$PlotName%in% grasslands,], zeros=T,envplot=envplot, var="SR", min.occur=10, min.class=3, alpha=0.01)
classSR.wood=classSR.fun(db=databp[databp$PlotName%in% woodlands,], zeros=T,envplot=envplot,var="SR",min.occur=10, min.class=3, alpha=0.01)
classSRnat.grass=classSR.fun(db=databp[databp$PlotName%in% grasslands,],  zeros=T,envplot=envplot,var="SRnat", min.occur=10, min.class=3, alpha=0.01)
classSRnat.wood=classSR.fun(db=databp[databp$PlotName%in% woodlands,], zeros=T,envplot=envplot,var="SRnat",min.occur=10, min.class=3, alpha=0.01)
classSRali.grass=classSR.fun(db=databp[databp$PlotName%in% grasslands,], zeros=T, envplot=envplot,var="SRali", min.occur=10, min.class=3, alpha=0.01)
classSRali.wood=classSR.fun(db=databp[databp$PlotName%in% woodlands,], zeros=T,envplot=envplot,var="SRali",min.occur=10, min.class=3, alpha=0.01)

#### summary of effects
effect.summary=data.frame(classSR, nat=classSRnat[,c('k.pval','rho', 'S.pval','effect')],
           ali=classSRali[,c('k.pval','rho', 'S.pval','effect')],
           wee=classSRwee[,c('k.pval','rho', 'S.pval','effect')])

effect.summary$total.effect=rowSums(abs(effect.summary[,c('effect','nat.effect','ali.effect')]))
                          

#######TABLE for supplementary material

tmp=effect.summary[, c("SpeciesCode","Division","Class","Family","SpeciesName","Growth.forms","ALIEN","WEEDOC",
               "Sp.occurence.x","nb.class","max","min",
               "effect","rho","S.pval",
               "nat.effect","nat.rho","nat.S.pval",
               "ali.effect","ali.rho","ali.S.pval",
                       "total.effect")]

write.csv(tmp, file="effect.summary.csv")


##################### STATS #################
# proportions of species effects for alien vs, native 
effect.table=as.data.frame(cbind(table(effect.summary$ALIEN,effect.summary$effect),
                          table(effect.summary$ALIEN,effect.summary$nat.effect),
                            table(effect.summary$ALIEN,effect.summary$ali.effect)))
names(effect.table)= c('Negative.SR','Neutral.SR','Positive.SR',
                       'Negative.SRnat','Neutral.SRnat','Positive.SRnat',
                       'Negative.SRali','Neutral.SRali','Positive.SRali')                      
rownames=c('Natives','Aliens')

#chisquare test on frequency of +, - and ns tests 
(Xsq <- chisq.test(effect.table[,1:3])) # barely significant
text(x=3.5, y=max(effect.table[,-c(2,5,8)]/sums*100)+0.5, label=(p2star(Xsq$p.val)))
(Xsq <- chisq.test(effect.table[,4:6])) # significant
text(x=9.5, y=max(effect.table[,-c(2,5,8)]/sums*100)+0.5, label=(p2star(Xsq$p.val)))
(Xsq <- chisq.test(effect.table[,7:9])) # significant
text(x=15.5, y=max(effect.table[,-c(2,5,8)]/sums*100)+0.5, label=(p2star(Xsq$p.val)))


# Differences between effect of alien on SR, SRnqt qnd SRali, as if 3 different samples :
tmp=matrix(effect.table[2,], 3,3, byrow=T)
mode(tmp)="numeric"
(Xsq <- chisq.test(tmp))   # significant difference of effects (on SRali)
Xsq$residuals

# Differences between effect of nqtives on SR, SRnat qnd SRali, as if 3 different samples :
tmp=matrix(effect.table[1,], 3,3, byrow=T)
mode(tmp)="numeric"
(Xsq <- chisq.test(tmp))   # significant difference of effects
Xsq$residuals


############### Overall distribution of rho values 
hist(effect.summary$rho, breaks=20)
abline(v=quantile(effect.summary$rho, 0.025, na.rm=T))
abline(v=quantile(effect.summary$rho, 0.975, na.rm=T))

hist(log(effect.summary$S.pval), breaks=40)
abline(v=log(quantile(effect.summary$S.pval, 0.025, na.rm=T)), col='red')
abline(v=log(quantile(effect.summary$S.pval, 0.975, na.rm=T)), col='red')

plot(effect.summary$rho,log(effect.summary$S.pval))
abline(v=quantile(effect.summary$rho, 0.025, na.rm=T))
abline(v=quantile(effect.summary$rho, 0.975, na.rm=T))
abline(h=log(quantile(effect.summary$S.pval, 0.025, na.rm=T)), col='red')
abline(h=log(0.01), col='red')

#### Comparing effects on AR and NR
tmp=effect.summary
#tmp=effect.summary[tmp$total.effect>0 & !is.na(tmp$total.effect),]
plot(-tmp$nat.rho,-tmp$ali.rho, ylim=c(-1,1), xlim=c(-1,1),
     col=c('forestgreen', 'firebrick')[as.numeric(tmp$ALIEN)+1],
     pch=c(21,20)[as.numeric(abs(tmp$total.effect)>0)+1])
summary(lm(nat.rho~ ali.rho*ALIEN, tmp))
summary(lm(nat.rho~ ali.rho, tmp))
abline(v=0, h=0)

tmp=grass.effect.summary
tmp=tmp[tmp$total.effect>0 & !is.na(tmp$total.effect),]
plot(-tmp$nat.rho,-tmp$ali.rho, ylim=c(-1,1), xlim=c(-1,1),
     col=c('forestgreen', 'firebrick')[as.numeric(tmp$ALIEN)+1],
     pch=c(21,20)[as.numeric(abs(tmp$total.effect)>0)+1])
summary(lm(nat.rho~ ali.rho*ALIEN, tmp))
summary(lm(nat.rho~ ali.rho, tmp))
abline(v=0, h=0)

tmp=wood.effect.summary
tmp=tmp[tmp$total.effect>0 & !is.na(tmp$total.effect),]
plot(-tmp$nat.rho,-tmp$ali.rho, ylim=c(-1,1), xlim=c(-1,1),
     col=c('forestgreen', 'firebrick')[as.numeric(tmp$ALIEN)+1],
     pch=c(21,20)[as.numeric(abs(tmp$total.effect)>0)+1])
summary(lm(nat.rho~ ali.rho*ALIEN, tmp))
summary(lm(nat.rho~ ali.rho, tmp))
abline(v=0, h=0)

plot(tmp$nat.effect,jitter(tmp$ali.effect,factor=0.2), ylim=c(-1.3,1.3), xlim=c(-1.2,1.2),pch=20,
     col=c('forestgreen', 'firebrick')[as.numeric(tmp$ALIEN)+1])

text(tmp$nat.effect,jitter(tmp$ali.effect, factor=2.2), 
     label=tmp$SpeciesCode, cex=0.7,
     col=c('forestgreen', 'firebrick')[as.numeric(tmp$ALIEN)+1])

table(tmp$ali.effect,tmp$nat.effect)

## PCA
rownames(effect.summary)=effect.summary$SpeciesCode
tmp[,c('nat.rho','rho','ali.rho')]=-tmp[,c('nat.rho','rho','ali.rho')]
aa=PCA(tmp[,c('ALIEN','nat.rho','rho','ali.rho')],quali.sup=1)

bb <- cbind.data.frame(tmp$ALIEN,aa$ind$coord)
ell <- coord.ellipse(bb,bary=TRUE,level.conf = 0.99)
     
plot(aa, choix="ind",habillage=1,col.hab=c('forestgreen', 'firebrick'),
     label="none", axes=c(1,2), xaxt='n', yaxt='n', ann=F)

# average effects :
t.test(-rho ~ ALIEN , classSR) # non significant
t.test(-rho ~ ALIEN , classSRali)   # signif
t.test(-rho ~ ALIEN , classSRnat)   # signif

tmp=boxplot(rho ~ ALIEN , classSR, varwidth=T, ylab='Spearman rho', xaxt='n', col=c('forestgreen', 'firebrick'))
tmp=boxplot(rho ~ ALIEN , classSRnat, varwidth=T, ylab='Spearman rho', xaxt='n', col=c('forestgreen', 'firebrick'))
tmp=boxplot(rho ~ ALIEN , classSRali, varwidth=T, ylab='Spearman rho', xaxt='n', col=c('forestgreen', 'firebrick'))

# Are all the rho values significantly positive or negative ??
t.test(classSR$rho[classSR$ALIEN==1],alternative = "greater") # not signif

# rho values BOXPLOT :
x11()
par(mfrow=c(1,3), mar=c(3,4,3,1), las=1, cex=0.9)

tmp=boxplot(-rho ~ ALIEN , classSR, varwidth=T, ylab="Spearman's rho",
            col=c('lightgrey', 'grey40'),ylim=c(-0.8,0.8),
            names=c('Natives', 'Aliens'))
abline(h=0)
title(main="a) Correlations with total richness", font=1, cex.main=0.8)
t=t.test(classSR$rho[classSR$ALIEN==1],alternative = "greater") # signif
text(x=2.2, y=0.2, label=(p2star(t$p.val)))
t=t.test(classSR$rho[classSR$ALIEN==0],alternative = "less") # signif
text(x=1.2, y=0.3, label=(p2star(t$p.val)))

tmp=boxplot(-rho ~ ALIEN , classSRnat, varwidth=T, ylab="Spearman's rho",
            col=c('lightgrey', 'grey40'),ylim=c(-0.8,0.8),
            names=c('Natives', 'Aliens'))
abline(h=0)
title(main="b) Correlations with Native richness", font=1, cex.main=0.8)
t=t.test(classSRnat$rho[classSRnat$ALIEN==1],alternative = "greater") # signif
text(x=2.2, y=0.2, label=(p2star(t$p.val)))
t=t.test(classSRnat$rho[classSRnat$ALIEN==0],alternative = "less") # signif
text(x=1.2, y=0.3, label=(p2star(t$p.val)))

tmp=boxplot(-rho ~ ALIEN , classSRali, varwidth=T,
            ylab="Spearman's rho", col=c('lightgrey', 'grey40'),ylim=c(-0.8,0.8),
            names=c('Natives', 'Aliens'))
abline(h=0)
title(main="c) Correlations with Alien richness", font=1,cex.main=0.8)
t=t.test(classSRali$rho[classSRali$ALIEN==0],alternative = "greater") # signif
text(x=1.2, y=0.2, label=(p2star(t$p.val)))
t=t.test(classSRali$rho[classSRali$ALIEN==1],alternative = "less") # signif
text(x=2.2, y=0.3, label=(p2star(t$p.val)))



# significant non zero rho values for effects on natives and aliens in GRASSLANDS:
x11()

par(mfrow=c(1,3), mar=c(3,4,3,1), las=1, cex=0.9)

tmp=boxplot(-rho ~ ALIEN , classSR.grass, varwidth=T, ylab="Spearman's rho",
            col=c('lightgrey', 'grey40'),ylim=c(-0.8,0.8),
            names=c('Natives', 'Aliens'))
abline(h=0)
title(main="a) Correlations with total richness", font=1, cex.main=0.8)
t=t.test(classSR.grass$rho[classSR.grass$ALIEN==1],alternative = "greater") # signif
text(x=2.2, y=0.2, label=(p2star(t$p.val)))
t=t.test(classSR.grass$rho[classSR.grass$ALIEN==0],alternative = "less") # signif
text(x=1.2, y=0.3, label=(p2star(t$p.val)))

tmp=boxplot(-rho ~ ALIEN , classSRnat.grass, varwidth=T, ylab="Spearman's rho",
            col=c('lightgrey', 'grey40'),ylim=c(-0.8,0.8),
            names=c('Natives', 'Aliens'))
abline(h=0)
title(main="b) Correlations with Native richness", font=1, cex.main=0.8)
t=t.test(classSRnat.grass$rho[classSRnat.grass$ALIEN==1],alternative = "greater") # signif
text(x=2.2, y=0.2, label=(p2star(t$p.val)))
t=t.test(classSRnat.grass$rho[classSRnat.grass$ALIEN==0],alternative = "less") # signif
text(x=1.2, y=0.3, label=(p2star(t$p.val)))

tmp=boxplot(-rho ~ ALIEN , classSRali.grass, varwidth=T,
            ylab="Spearman's rho", col=c('lightgrey', 'grey40'),ylim=c(-0.8,0.8),
            names=c('Natives', 'Aliens'))
abline(h=0)
title(main="c) Correlations with Alien richness", font=1,cex.main=0.8)
t=t.test(classSRali.grass$rho[classSRali.grass$ALIEN==0],alternative = "greater") # signif
text(x=1.2, y=0.2, label=(p2star(t$p.val)))
t=t.test(classSRali.grass$rho[classSRali.grass$ALIEN==1],alternative = "less") # signif
text(x=2.2, y=0.3, label=(p2star(t$p.val)))

# significant non zero rho values for effects on natives and aliens in woodLANDS:
x11()

par(mfrow=c(1,3), mar=c(3,4,3,1), las=1, cex=0.9)

tmp=boxplot(-rho ~ ALIEN , classSR.wood, varwidth=T, ylab="Spearman's rho",
            col=c('lightgrey', 'grey40'),ylim=c(-0.8,0.8),
            names=c('Natives', 'Aliens'))
abline(h=0)
title(main="a) Correlations with total richness", font=1, cex.main=0.8)
t=t.test(classSR.wood$rho[classSR.wood$ALIEN==1],alternative = "greater") # signif
text(x=2.2, y=0.2, label=(p2star(t$p.val)))
t=t.test(classSR.wood$rho[classSR.wood$ALIEN==0],alternative = "less") # signif
text(x=1.2, y=0.3, label=(p2star(t$p.val)))

tmp=boxplot(-rho ~ ALIEN , classSRnat.wood, varwidth=T, ylab="Spearman's rho",
            col=c('lightgrey', 'grey40'),ylim=c(-0.8,0.8),
            names=c('Natives', 'Aliens'))
abline(h=0)
title(main="b) Correlations with Native richness", font=1, cex.main=0.8)
t=t.test(classSRnat.wood$rho[classSRnat.wood$ALIEN==1],alternative = "greater") # signif
text(x=2.2, y=0.2, label=(p2star(t$p.val)))
t=t.test(classSRnat.wood$rho[classSRnat.wood$ALIEN==0],alternative = "less") # signif
text(x=1.2, y=0.3, label=(p2star(t$p.val)))

tmp=boxplot(-rho ~ ALIEN , classSRali.wood, varwidth=T,
            ylab="Spearman's rho", col=c('lightgrey', 'grey40'),ylim=c(-0.8,0.8),
            names=c('Natives', 'Aliens'))
abline(h=0)
title(main="c) Correlations with Alien richness", font=1,cex.main=0.8)
t=t.test(classSRali.wood$rho[classSRali.wood$ALIEN==0],alternative = "greater") # signif
text(x=1.2, y=0.2, label=(p2star(t$p.val)))
t=t.test(classSRali.wood$rho[classSRali.wood$ALIEN==1],alternative = "less") # signif
text(x=2.2, y=0.3, label=(p2star(t$p.val)))


# plot distribution of SR for each species dom class
# effect of aliens on SR
tmp=databp[databp$Sp.occurence>10,]
tmp=tmp[tmp$SpeciesCode %in% effect.summary$SpeciesCode[effect.summary$total.effect>0 & effect.summary$ALIEN==1],]
x11()
par(mfrow=c(4,4), mar=c(3,2,2,1))
for (i in unique(tmp$SpeciesCode) ) {
  x=tmp[as.character(tmp$SpeciesCode)==i,]
  if(length(unique(x$domlevels))<2) plot(SR~as.factor(domlevels), x, border="grey")
  if(length(unique(x$domlevels))>2){
  co=cor.test(x$SR,x$domlevels, method='spearman')
  #co=cor.test(x$SR,x$domlevels, method='kendall')
 
#   k=kruskal.test(x$SR,x$domlevels)
#   
#   if (k$p.val<0.05 & co$est<0) plot(SR~as.factor(domlevels), x, col="forestgreen")
#   if (k$p.val<0.05 & co$est>0) plot(SR~as.factor(domlevels), x, col="firebrick")
#   if (k$p.val>0.05) plot(SR~as.factor(domlevels), x, col="grey")
  
  if (co$p.val<0.05 & co$est<0) plot(SR~as.factor(domlevels), x, col="forestgreen")
  if (co$p.val<0.05 & co$est>0) plot(SR~as.factor(domlevels), x, col="firebrick")
  if (co$p.val>0.05) plot(SR~as.factor(domlevels), x, col="grey")

  title(main=as.character(i))
  mtext(3,text=round(co$est,2), adj=1,cex=0.6)
  }
  }
  
# correlation between effect on natives and on aliens :
plot(nat.rho~ali.rho, effect.summary,pch=20,
col=c('forestgreen', 'firebrick')[as.numeric(effect.summary$ALIEN==1)+1],
     xlab='effect on aliens (rho)', ylab='effect on natives') 
fit=lm(nat.rho~ali.rho + ALIEN, effect.summary)
summary(fit)
abline(fit)
abline(h=0, col='grey')
abline(v=0, col='grey')

plot(nat.rho~ali.rho, wood.effect.summary,pch=20,
col=c('forestgreen', 'firebrick')[as.numeric(effect.summary$ALIEN==1)+1],
     xlab='effect on aliens (rho)', ylab='effect on natives') 
fit=lm(nat.rho~ali.rho + ALIEN, effect.summary)
summary(fit)
abline(fit)
abline(h=0, col='grey')
abline(v=0, col='grey')

# 
# 
# ############ Effect of Life duration ######
# 
# hist(effect.summary$rho, breaks=15)
# hist(effect.summary$rho[effect.summary$LIFE_DURATION=='P'], col='firebrick', breaks=15, add=T)
# hist(effect.summary$rho[effect.summary$LIFE_DURATION=='A'], col='forestgreen', breaks=15, add=T)
# 
# # boxplot of rho values :
# 
# tmp=boxplot(-rho ~ as.character(lifestyle) , classSR, varwidth=T, xaxt='n',lab='-Spearman rho', col=c('forestgreen', 'firebrick'))
# 
# tmp=boxplot(-rho ~ as.character(lifestyle) + effect , classSR, varwidth=T, xaxt='n',lab='-Spearman rho', col=c('forestgreen', 'firebrick'))
# axis(1, at=1:6, labels = rep(c('A','P'),3), tcl=-0.3)
# axis(1, at=c(1.5,3.5,5.5), line=1.5,lty=0,labels = c('Negative effect','Neutral', 'Positive effect'))
# axis(1, at=c(0.5,2.5,4.5,6.5), tcl=-4,labels = rep('',4))
# mtext(side=1,at=1:6,text=tmp$n, cex=0.75, line=0.2)
# abline(h=0)
# 
# tmp=boxplot(-rho ~ as.character(lifestyle) , classSRnat, varwidth=T, xaxt='n',lab='-Spearman rho', col=c('forestgreen', 'firebrick'))
# 
# tmp=boxplot(-rho ~ as.character(lifestyle) + effect , classSRnat, varwidth=T, xaxt='n',lab='-Spearman rho', col=c('forestgreen', 'firebrick'))
# axis(1, at=1:6, labels = rep(c('A','P'),3), tcl=-0.3)
# axis(1, at=c(1.5,3.5,5.5), line=1.5,lty=0,labels = c('Negative effect','Neutral', 'Positive effect'))
# axis(1, at=c(0.5,2.5,4.5,6.5), tcl=-4,labels = rep('',4))
# mtext(side=1,at=1:6,text=tmp$n, cex=0.75, line=0.2)
# abline(h=0)
# 
# 
# # proportions of species effects for annuals vs. perennial
# effect.table=as.data.frame(cbind(table(effect.summary$lifestyle,effect.summary$effect),
#                           table(effect.summary$lifestyle,effect.summary$nat.effect),
#                             table(effect.summary$lifestyle,effect.summary$ali.effect)))
# names(effect.table)= c('Negative.SR','Neutral.SR','Positive.SR',
#                        'Negative.SRnat','Neutral.SRnat','Positive.SRnat',
#                        'Negative.SRali','Neutral.SRali','Positive.SRali')                      
# rownames=c('Annuals','Perennial')
# 
# #chisquare test on frequency of +, - and ns tests
# summary(table(lifestyle=effect.summary$lifestyle,effect=effect.summary$effect,type=effect.summary$ALIEN))
# summary(table(lifestyle=effect.summary$lifestyle,effect=effect.summary$effect))
# summary(table(type=effect.summary$ALIEN,effect=effect.summary$effect))
#                        
# summary(table(lifestyle=effect.summary$lifestyle,effect=effect.summary$nat.effect,type=effect.summary$ALIEN))
# summary(table(lifestyle=effect.summary$lifestyle,effect=effect.summary$nat.effect))
# 
# summary(table(lifestyle=effect.summary$lifestyle,effect=effect.summary$ali.effect,type=effect.summary$ALIEN))
# summary(table(lifestyle=effect.summary$lifestyle,effect=effect.summary$ali.effect))
# 
# (Xsq <- chisq.test(effect.table[,1:3])) # not significant
# (Xsq <- chisq.test(effect.table[,4:6])) # not significant
# (Xsq <- chisq.test(effect.table[,7:9])) # not significant
# 
# 
# barplot(as.matrix(effect.table[,-c(2,5,8)]/c(58,144)*100),
#         beside=T,names=F, cex.names=0.8, col=c('green', 'darkgreen'))
# axis(1, at=c(2,5,8,11,14,17), line=-0.5, labels = rep(c('-','+'),3), tcl=-0.3, lty=0, cex.axis=2)
# axis(1, at=c(3.5,9.5,15.5), line=0.5,lty=0,labels = c('Total richness','Native richness', 'Alien richness'))
# axis(2, las=1)
# box(bty='l')
# legend(x='top', legend=c('Annual','Perennial'), fill=c('green', 'darkgreen'), cex=0.8)                        
# title(ylab='% species')
                       

##### proportions of species effects for growthforms
effect.table=as.data.frame(cbind(table(effect.summary$Growth.forms,effect.summary$effect),
                          table(effect.summary$Growth.forms,effect.summary$nat.effect),
                            table(effect.summary$Growth.forms,effect.summary$ali.effect)))
names(effect.table)= c('Negative.SR','Neutral.SR','Positive.SR',
                       'Negative.SRnat','Neutral.SRnat','Positive.SRnat',
                       'Negative.SRali','Neutral.SRali','Positive.SRali')                      
rownames=c('FE','GR','HR','SH','TR')
#fern, grass,herb,shru,tree

# Effect of grozth forms alone
(Xsq <- chisq.test(effect.table[,1:3])) # not significant
(Xsq <- chisq.test(effect.table[,4:6])) # significant
(Xsq <- chisq.test(effect.table[,7:9])) # significant

# Barplots of proportions :
par(mfrow=c(3,1), mar=c(3,4,1,1))
tmp=effect.table[,1:3]
tmp=tmp/rowSums( tmp)*100
barplot(as.matrix(tmp),
        beside=T, cex.names=1, col=colgrf)
box(bty='l')
legend(x='topright', legend=c('FE','GR','HR','SH','TR'), fill=colgrf, cex=0.9)                        
title(ylab='% species')

tmp=effect.table[,4:6]
tmp=tmp/rowSums( tmp)*100
barplot(as.matrix(tmp),
        beside=T, cex.names=1, col=colgrf)
box(bty='l')
title(ylab='% species')

tmp=effect.table[,7:9]
tmp=tmp/rowSums( tmp)*100
barplot(as.matrix(tmp),
        beside=T, cex.names=1, col=colgrf)
box(bty='l')
title(ylab='% species')



# Interaction between alien status and grf 

# EFFECt on SR
summary(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$ALIEN,effect=effect.summary$effect))
ftable(table(type=effect.summary$ALIEN,lifestyle=effect.summary$Growth.forms,effect=effect.summary$effect))

summary(table(lifestyle=effect.summary$Growth.forms,effect=effect.summary$effect))
ftable(table(lifestyle=effect.summary$Growth.forms,effect=effect.summary$effect))
summary(table(type=effect.summary$ALIEN,effect=effect.summary$effect))

#Test more specifically within negative effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==-1],type=effect.summary$ALIEN[effect.summary$effect==-1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==-1],type=effect.summary$ALIEN[effect.summary$effect==-1]), beside=T)
#Test more specifically within positive effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==1],type=effect.summary$ALIEN[effect.summary$effect==1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==1],type=effect.summary$ALIEN[effect.summary$effect==1]), beside=T)
#Test more specifically within neutral effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==0],type=effect.summary$ALIEN[effect.summary$effect==0]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==0],type=effect.summary$ALIEN[effect.summary$effect==0]), beside=T)


x11()
par(mar=c(3,5,1,1))
barplot(ftable(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$ALIEN,effect=effect.summary$effect)),
        beside=T, density=c(30,-1), offset=0, horiz=T,
        col=as.vector(rbind(colgrf, colgrf)))

legend(50,30, legend=c('Fern','Grass','Herb','Shrub','Tree'),
       fill=c("brown","goldenrod","palegreen","forestgreen","darkgreen"),
       bty='n',cex=0.8)
legend(35,28, legend=c("native","alien"),fill='darkgrey',density=c(30,-1), bty='n',cex=0.8)
axis(2, at=c(6,17,30), line=1.3,lty=0, adj=0.5, labels = c('-','neutral', '+'))
axis(2, at=c(rep(seq(2,10,2),3)+c(rep(0,5),rep(11,5),rep(22,5))) , 
     line=-1,lty=0, adj=0.5, las=2,cex.axis=0.8,
     labels = rep(c('Fern','Grass','Herb','Shrub','Tree'),3))
     axis(2, at=c(rep(seq(1,11,2),3)+c(rep(0,6),rep(11,6),rep(22,6))) , 
     lty=1,tcl=-0.2,line=-0.25,lwd.ticks=1, lwd=0,
     labels =rep('',18))
title(ylab='Effect on total richness', line=4)

          
# EFFECt on SR natives
summary(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$ALIEN,effect=effect.summary$nat.effect))
ftable(table(type=effect.summary$ALIEN,lifestyle=effect.summary$Growth.forms,effect=effect.summary$nat.effect))

summary(table(lifestyle=effect.summary$Growth.forms,effect=effect.summary$nat.effect))
summary(table(type=effect.summary$ALIEN,effect=effect.summary$nat.effect))

#Test more specifically within negative effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==-1],type=effect.summary$ALIEN[effect.summary$nat.effect==-1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==-1],type=effect.summary$ALIEN[effect.summary$nat.effect==-1]), beside=T)
#Test more specifically within positive effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==1],type=effect.summary$ALIEN[effect.summary$nat.effect==1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==1],type=effect.summary$ALIEN[effect.summary$nat.effect==1]), beside=T)
#Test more specifically within neutral effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==0],type=effect.summary$ALIEN[effect.summary$nat.effect==0]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==0],type=effect.summary$ALIEN[effect.summary$nat.effect==0]), beside=T)


x11()
par(mar=c(3,5,1,1))
barplot(ftable(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$ALIEN,effect=effect.summary$nat.effect)),
        beside=T, density=c(30,-1), offset=0, horiz=T,
        col=as.vector(rbind(colgrf, colgrf)))

legend(50,30, legend=c('Fern','Grass','Herb','Shrub','Tree'),
       fill=c("brown","goldenrod","palegreen","forestgreen","darkgreen"),
       bty='n',cex=0.8)
legend(35,28, legend=c("native","weed"),fill='darkgrey',density=c(30,-1), bty='n',cex=0.8)
axis(2, at=c(6,17,30), line=1.3,lty=0, adj=0.5, labels = c('-','neutral', '+'))
axis(2, at=c(rep(seq(2,10,2),3)+c(rep(0,5),rep(11,5),rep(22,5))) , 
     line=-1,lty=0, adj=0.5, las=2,cex.axis=0.8,
     labels = rep(c('Fern','Grass','Herb','Shrub','Tree'),3))
axis(2, at=c(rep(seq(1,11,2),3)+c(rep(0,6),rep(11,6),rep(22,6))) , 
     lty=1,tcl=-0.2,line=-0.25,lwd.ticks=1, lwd=0,
     labels =rep('',18))
title(ylab='Effect on Native richness', line=4)

# EFFECt on SR aliens
summary(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$ALIEN,effect=effect.summary$ali.effect))
ftable(table(type=effect.summary$ALIEN,lifestyle=effect.summary$Growth.forms,effect=effect.summary$ali.effect))

summary(table(lifestyle=effect.summary$Growth.forms,effect=effect.summary$ali.effect))
summary(table(type=effect.summary$ALIEN,effect=effect.summary$ali.effect))

#Test more specifically within negative effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$ali.effect==-1],type=effect.summary$ALIEN[effect.summary$ali.effect==-1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$ali.effect==-1],type=effect.summary$ALIEN[effect.summary$ali.effect==-1]), beside=T)
#Test more specifically within positive effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$ali.effect==1],type=effect.summary$ALIEN[effect.summary$ali.effect==1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$ali.effect==1],type=effect.summary$ALIEN[effect.summary$ali.effect==1]), beside=T)
#Test more specifically within neutral effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$ali.effect==0],type=effect.summary$ALIEN[effect.summary$ali.effect==0]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$ali.effect==0],type=effect.summary$ALIEN[effect.summary$ali.effect==0]), beside=T)


x11()
par(mar=c(3,5,1,1))
barplot(ftable(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$ALIEN,effect=effect.summary$ali.effect)),
        beside=T, density=c(30,-1), offset=0, horiz=T,
        col=as.vector(rbind(colgrf, colgrf)))

legend(50,30, legend=c('Fern','Grass','Herb','Shrub','Tree'),
       fill=c("brown","goldenrod","palegreen","forestgreen","darkgreen"),
       bty='n',cex=0.8)
legend(35,28, legend=c("native","alien"),fill='darkgrey',density=c(30,-1), bty='n',cex=0.8)
axis(2, at=c(6,17,30), line=1.3,lty=0, adj=0.5, labels = c('-','neutral', '+'))
axis(2, at=c(rep(seq(2,10,2),3)+c(rep(0,5),rep(11,5),rep(22,5))) , 
     line=-1,lty=0, adj=0.5, las=2,cex.axis=0.8,
     labels = rep(c('Fern','Grass','Herb','Shrub','Tree'),3))
axis(2, at=c(rep(seq(1,11,2),3)+c(rep(0,6),rep(11,6),rep(22,6))) , 
     lty=1,tcl=-0.2,line=-0.25,lwd.ticks=1, lwd=0,
     labels =rep('',18))
title(ylab='Effect on Alien richness', line=4)


############  Interaction between WEEDOC status and grf 
table(alien=effect.summary$ALIEN,weedoc=effect.summary$WEEDOC)

# EFFECt on SR
summary(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$WEEDOC,effect=effect.summary$effect))
summary(table(lifestyle=effect.summary$Growth.forms,effect=effect.summary$effect))
summary(table(type=effect.summary$WEEDOC,effect=effect.summary$effect))

#Test more specifically within negative effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==-1],type=effect.summary$WEEDOC[effect.summary$effect==-1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==-1],type=effect.summary$WEEDOC[effect.summary$effect==-1]), beside=T)
#Test more specifically within positive effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==1],type=effect.summary$WEEDOC[effect.summary$effect==1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==1],type=effect.summary$WEEDOC[effect.summary$effect==1]), beside=T)
#Test more specifically within neutral effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==0],type=effect.summary$WEEDOC[effect.summary$effect==0]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$effect==0],type=effect.summary$WEEDOC[effect.summary$effect==0]), beside=T)


x11()
par(mar=c(3,5,1,1))
barplot(ftable(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$WEEDOC,effect=effect.summary$effect)),
        beside=T, density=c(30,-1), offset=0, horiz=T,
        col=as.vector(rbind(colgrf, colgrf)))

legend(50,30, legend=c('Fern','Grass','Herb','Shrub','Tree'),
       fill=c("brown","goldenrod","palegreen","forestgreen","darkgreen"),
       bty='n',cex=0.8)
legend(35,28, legend=c("native","weed"),fill='darkgrey',density=c(30,-1), bty='n',cex=0.8)
axis(2, at=c(6,17,30), line=1.3,lty=0, adj=0.5, labels = c('-','neutral', '+'))
axis(2, at=c(rep(seq(2,10,2),3)+c(rep(0,5),rep(11,5),rep(22,5))) , 
     line=-1,lty=0, adj=0.5, las=2,cex.axis=0.8,
     labels = rep(c('Fern','Grass','Herb','Shrub','Tree'),3))
     axis(2, at=c(rep(seq(1,11,2),3)+c(rep(0,6),rep(11,6),rep(22,6))) , 
     lty=1,tcl=-0.2,line=-0.25,lwd.ticks=1, lwd=0,
     labels =rep('',18))
title(ylab='Effect on total richness', line=4)


# EFFECt on SR natives
summary(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$WEEDOC,effect=effect.summary$nat.effect))
summary(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$WEEDOC,effect=effect.summary$nat.effect))

summary(table(lifestyle=effect.summary$Growth.forms,effect=effect.summary$nat.effect))
summary(table(type=effect.summary$WEEDOC,effect=effect.summary$nat.effect))

#Test more specifically within negative effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==-1],type=effect.summary$WEEDOC[effect.summary$nat.effect==-1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==-1],type=effect.summary$WEEDOC[effect.summary$nat.effect==-1]), beside=T)
#Test more specifically within positive effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==1],type=effect.summary$WEEDOC[effect.summary$nat.effect==1]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==1],type=effect.summary$WEEDOC[effect.summary$nat.effect==1]), beside=T)
#Test more specifically within neutral effects
summary(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==0],type=effect.summary$WEEDOC[effect.summary$nat.effect==0]))
barplot(table(lifestyle=effect.summary$Growth.forms[effect.summary$nat.effect==0],type=effect.summary$WEEDOC[effect.summary$nat.effect==0]), beside=T)


x11()
par(mar=c(3,5,1,1))
barplot(ftable(table(lifestyle=effect.summary$Growth.forms,type=effect.summary$WEEDOC,effect=effect.summary$nat.effect)),
        beside=T, density=c(30,-1), offset=0, horiz=T,
        col=as.vector(rbind(colgrf, colgrf)))

legend(50,30, legend=c('Fern','Grass','Herb','Shrub','Tree'),
       fill=c("brown","goldenrod","palegreen","forestgreen","darkgreen"),
       bty='n',cex=0.8)
legend(30,28, legend=c("non-weed","weed"),fill='darkgrey',density=c(30,-1), bty='n',cex=0.8)
axis(2, at=c(6,17,30), line=1.3,lty=0, adj=0.5, labels = c('-','neutral', '+'))
axis(2, at=c(rep(seq(2,10,2),3)+c(rep(0,5),rep(11,5),rep(22,5))) , 
     line=-1,lty=0, adj=0.5, las=2,cex.axis=0.8,
     labels = rep(c('Fern','Grass','Herb','Shrub','Tree'),3))
axis(2, at=c(rep(seq(1,11,2),3)+c(rep(0,6),rep(11,6),rep(22,6))) , 
     lty=1,tcl=-0.2,line=-0.25,lwd.ticks=1, lwd=0,
     labels =rep('',18))
title(ylab='Effect on Native richness', line=4)


