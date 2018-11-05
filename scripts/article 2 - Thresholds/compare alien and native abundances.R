# Compare abundance distribtuions of native and aliens
db <- databp[databp$PlotName %in% unimprovedgrasslands,]


par(mfrow = c(1,2))
barplot(table(db$abun[which(db$NATIVE ==TRUE)]), main = "Natives", xlab = "abundance classes")
barplot(table(db$abun[which(db$ALIEN ==TRUE)]), main = "Aliens", xlab = "abundance classes" )

barplot(table(db$abun[which(db$NATIVE ==TRUE)])/sum(db$NATIVE ==TRUE), main = "Natives", xlab = "abundance classes", ylim = c(0,0.60), ylab = "Relative frequency")
barplot(table(db$abun[which(db$ALIEN ==TRUE)])/sum(db$ALIEN ==TRUE, na.rm = T), main = "Aliens", xlab = "abundance classes", ylim = c(0,0.60))

tmp  = t(as.matrix(table(db$abun, db$NATIVE)))/as.vector(table(db$NATIVE))
chisq.test(t(tmp))

## test within subsets with target species

par(mfrow = c(2,4), xpd = T, mar = c(3,3,1,1))
for (sp in impsp){
  dat <- db[ db$SpeciesCode == sp,]
  ab <- unique(dat$abun)

  plot( c(0,7),c(0,7), type = "n") 
  for (i in ab) {
    subs <- db[ db$PlotName %in% unique( dat[dat$abun == i,]$PlotName),]
    subs <- subs[subs$NATIVE == 1,]
    subs$abun <- factor(subs$abun, levels = 1:7)
    # boxplot( x =as.numeric(subs$abun), add= T, at = i)
    points( x = rep(i,7), y = 1:7, pch = 20, cex = (table(subs$abun)/length(subs$abun))*10, col = "firebrick")
    mtext(3,text = sp, cex = 0.7)
  }
}


par(mfrow = c(2,4), xpd =FALSE, mar = c(3,3,1,1))
for (sp in impsp){
  dat <- db[ db$SpeciesCode == sp,]
  subs <- db[ db$PlotName %in% unique(dat$PlotName),]
  subs <- subs[subs$NATIVE == 1,]
  barplot(table(subs$abun)/sum(subs$NATIVE), ylim = c(0,0.60), main = sp)
 abline(v = table.div.part[rownames(table.div.part) == sp,"th.CI"], col = "red")
}



### Barplots of mean distribution of abundances within plots
par(mfrow = c(2,4), xpd =FALSE, mar = c(3,3,1,1))
  subs <- comm[  ,names(comm) %in% natives]
  
  tbl <- t(apply(subs, 1, FUN = function(j) {
    return(table(factor(x= j, levels = 1:7)))
    }))
  
  for (sp in impsp){
  abun <- comm[ ,sp]
  m <- t(sapply(1:7, FUN = function(i) {
    return(if(sum(abun == i, na.rm = T)>0) (apply(tbl[abun == i,], 2, mean, na.rm = T)) else rep(NA,7))
  }))
  # serr <-t(sapply(1:7, FUN = function(i) apply(tbl[abun == i,], 2, sd)/nrow(tbl[abun == i,])))
  sdev <- t(sapply(1:7, FUN = function(i) {
    return(if(sum(abun == i, na.rm = T)>0) (apply(tbl[abun == i,], 2, sd, na.rm = T)) else rep(NA,7))
  }))
  b <- barplot(m, beside = T, ylim = c(0,max(m, na.rm = T) + max(sdev, na.rm = T))  )
  segments( b,  m +sdev,  b, m )
  
  abline(v = b[table.div.part[rownames(table.div.part) == sp,"th.CI"],table.div.part[rownames(table.div.part) == sp,"th.CI"]], col = "red")
  }
  
  
  
  ### Barplots of mean distribution of abundances within plots
  par(mfrow = c(1,2), xpd =FALSE, mar = c(3,3,1,1))
  subs <- comm[  ,names(comm) %in% natives]
  tbl <- t(apply(subs, 1, FUN = function(j) {
    return(table(factor(x= j, levels = 1:7)))
  }))
  
    m <- apply(tbl, 2, mean, na.rm = T)
    sdev <- apply(tbl, 2, sd, na.rm = T)
    b <- barplot(m, beside = T, ylim = c(0,10)  )
    segments( b,  m +sdev,  b, m )
    
    subs <- comm[  ,names(comm) %in% aliens]
    tbl <- t(apply(subs, 1, FUN = function(j) {
      return(table(factor(x= j, levels = 1:7)))
    }))
    
    m <- apply(tbl, 2, mean, na.rm = T)
    sdev <- apply(tbl, 2, sd, na.rm = T)
b <- barplot(m, beside = T, ylim =  c(0,10) , col = "firebrick"  )
segments( b,  m +sdev,  b, m )

