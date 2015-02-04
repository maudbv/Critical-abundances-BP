# Species richness - species dominance class relationships

################## ALL PLOTS ############################################


## Testing the relationships for each species with kruskal wallis + spearman
classSR=classSR.fun(db=databp, var='SR',zeros=F,env=envplot, min.occur=10, min.class=3, alpha=0.01)
classSRnat=classSR.fun(db=databp, zeros=F, var='SRnat',min.occur=10, env=envplot,min.class=3, alpha=0.01)
classSRali=classSR.fun(db=databp, zeros=F, var='SRali', min.occur=10, env=envplot,min.class=3, alpha=0.01)
classSRwee=classSR.fun(db=databp, zeros=F, var='SRwee', min.occur=10,env=envplot, min.class=3, alpha=0.01)

#### summary of effects
effect.summary=data.frame(classSR, nat=classSRnat[,c('k.pval','rho', 'S.pval','effect')],
                          ali=classSRali[,c('k.pval','rho', 'S.pval','effect')],
                          wee=classSRwee[,c('k.pval','rho', 'S.pval','effect')])

effect.summary$total.effect=rowSums(abs(effect.summary[,c('effect','nat.effect','ali.effect')]))

rownames(effect.summary) <- as.character(effect.summary$SpeciesCode )

#######TABLE for supplementary material

tmp=effect.summary[, c("SpeciesCode","Division","Class","Family","SpeciesName","Growth.forms","ALIEN","WEEDOC",
                       "Sp.occurence.x","nb.class","max","min",
                       "effect","rho","S.pval",
                       "nat.effect","nat.rho","nat.S.pval",
                       "ali.effect","ali.rho","ali.S.pval",
                       "total.effect")]

write.csv(tmp, file="effect.summary.csv")

##################### STATS across allplots #################
# proportions of species effects for alien vs, native 
effect.table=as.data.frame(cbind(table(effect.summary$ALIEN,effect.summary$effect),
                                 table(effect.summary$ALIEN,effect.summary$nat.effect),
                                 table(effect.summary$ALIEN,effect.summary$ali.effect)))
names(effect.table)= c('Negative.SR','Neutral.SR','Positive.SR',
                       'Negative.SRnat','Neutral.SRnat','Positive.SRnat',
                       'Negative.SRali','Neutral.SRali','Positive.SRali')                      
rownames=c('Natives','Aliens')


################## GRASSLANDS ##########################################
db=databp[databp$PlotName%in% grasslands,]

## Testing the relationships for each species with kruskal wallis + spearman
grass.classSR=classSR.fun(db=db, var='SR',zeros=F,  min.occur=10, min.class=3, alpha=0.01)
grass.classSRnat=classSR.fun(db=db, var='SRnat',zeros=F, min.occur=10, min.class=3, alpha=0.01)
grass.classSRali=classSR.fun(db=db, var='SRali', zeros=F, min.occur=10, min.class=3, alpha=0.01)

#### summary of effects
grass.effect.summary=data.frame(grass.classSR, nat=grass.classSRnat[,c('k.pval','rho', 'S.pval','effect')],
           ali=grass.classSRali[,c('k.pval','rho', 'S.pval','effect')])
grass.effect.summary$total.effect=rowSums(abs(grass.effect.summary[,c('effect','nat.effect','ali.effect')]))
grass.effect.summary=na.omit(grass.effect.summary)       
rownames(grass.effect.summary) <- as.character(grass.effect.summary$SpeciesCode )

#######TABLE for supplementary material
tmp=grass.effect.summary[, c("SpeciesCode","Division","Class","Family","SpeciesName","Growth.forms","ALIEN","WEEDOC",
               "Sp.occurence.x","nb.class","max","min",
               "effect","rho","S.pval",
               "nat.effect","nat.rho","nat.S.pval",
               "ali.effect","ali.rho","ali.S.pval",
                       "total.effect")]

write.csv(tmp, file="grass.effect.summary.csv")


################################ WOODLANDS

db=databp[databp$PlotName%in% woodlands,]

## Testing the relationships for each species with kruskal wallis + spearman
wood.classSR=classSR.fun(db=db, var='SR', min.occur=10, min.class=3, alpha=0.01)
wood.classSRnat=classSR.fun(db=db, var='SRnat',  min.occur=10, min.class=3, alpha=0.01)
wood.classSRali=classSR.fun(db=db, var='SRali', min.occur=10, min.class=3, alpha=0.01)

#### summary of effects
wood.effect.summary=data.frame(wood.classSR, nat=wood.classSRnat[,c('k.pval','rho', 'S.pval','effect')],
                               ali=wood.classSRali[,c('k.pval','rho', 'S.pval','effect')])

wood.effect.summary$total.effect=rowSums(abs(wood.effect.summary[,c('effect','nat.effect','ali.effect')]))
wood.effect.summary=na.omit(wood.effect.summary)       
rownames(wood.effect.summary) <- as.character(wood.effect.summary$SpeciesCode )


#######TABLE for supplementary material
tmp=wood.effect.summary[, c("SpeciesCode","Division","Class","Family","SpeciesName","Growth.forms","ALIEN","WEEDOC",
                            "Sp.occurence.x","nb.class","max","min",
                            "effect","rho","S.pval",
                            "nat.effect","nat.rho","nat.S.pval",
                            "ali.effect","ali.rho","ali.S.pval",
                            "total.effect")]
write.csv(tmp, file="wood.effect.summary.csv")


rm(db)