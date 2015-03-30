<<<<<<< HEAD
### Threshold analysis plotting functions

## plot summary outputs
plot.effect.summary <- function(  effects =effects.glmSR.grass ) {       
  
  n <- length(unique(effects$group)) 
  
  par(mfrow=c(n,2), oma=c(2,5,2,0), mar=c(2,2,1,1))
  
  for (i in 1:n) {
    S <- as.data.frame(effects[effects$group == unique(effects$group)[i],])
    barplot(S$prop.impact*100, ylim=c(0,50))
    if (i==1) mtext(side=3, text="% negative\neffects",adj=0.5,line=1, cex=0.8)
    
    barplot(S$prop.thr*100, col="black" , ylim=c(0,50))
    if (i==1) mtext(side=3, text="% threshold\neffects", adj=0.5,line=1, cex=0.8)
    
  }   
  mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=2.5, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
  mtext(text="Abundance class", side=1, outer=T, line=0.5)
}

### plotting individual species
plot.sp.glm=function(sp="ACHMIL", M=glmSR, var="SR", db=databp, type="boot") {  
  Z= db[db$SpeciesCode==sp,]
  if (type == "boot") M$thresh <-M$boot.thresh
  col=rep("white",6)
  col[ as.numeric( M$thresh[sp,"th"])]="red"
  boxplot(as.formula(paste(var, " ~ abun")), data=Z, xlim=c(0,6), outline=F, varwidth=T, col=col)
  mtext(side=3, line=-1.1, at=2:6,text=sapply(M$P[sp,], FUN=p2star), cex=0.8)
  mtext(side=3, text=paste ("rho" ,round(M$spearman[sp,"rho"],2), p2star(M$spearman[sp,"p.val"])), adj=1, cex=0.7)
  title(main=sp, cex.main=0.8, line=1)
  
  return(M$thresh[sp,"th"])
}

## plotting all significant species
plot.glm=function(M=glmSR, var="SR", db=databp, sel.criteria = c("spear")) {
  
  par(mfrow=c(5,5), mar=c(2,2,2,2))
  
  # select only species according to seletion criteria :
  if (sel.criteria == "spear") sel <- rownames(M$spearman) [M$spearman[,"p.val"] <=0.05 &  M$spearman[,"rho"] <0 ]
  if (sel.criteria == "dev1") sel <-  rownames(M$glms) [M$glms[,"dev.ratio"] > 0.05 &  !is.na(M$thresh[,"th"]) ]
  if (sel.criteria == "none") sel <- TRUE
  if (sel.criteria == "th.exist") sel <-  rownames(M$thresh) [ !is.na(M$thresh[,"th"]) ]
  
  ### Loop on selected species 
  for (sp in sel)  {
    th=plot.sp.glm(sp= sp, M=M, var=var, db=db)
    print(paste(sp,":",th))
  }
}

### thredhold barplots
thresh.prop=function(effects=effects.glmSR.grass, data=species, y=T,ylim=c(0,0.35)) {
  #graphical representation
  barplot(effects$prop.thr[effects$group=="ALIEN:0"], names.arg= paste("c",2:6, sep=""), 
          col="grey", ylim=ylim, las=1 )
  if (y==T) mtext(side=2, text="Prop. threshold effects",line=3.4, cex=0.8)
  barplot( effects$prop.thr[effects$group=="ALIEN:1"], names.arg= paste("c",2:6, sep=""),
           col="black" , ylim=ylim , las=1)
  if (y==T) mtext(side=2, text="Prop. threshold effects", line=3.4, cex=0.8)
  
  print(wilcox.test(effects$prop.thr[effects$group=="ALIEN:1"],effects$prop.thr[effects$group=="ALIEN:0"], paired=T))
  print(friedman.test(cbind(effects$prop.thr[effects$group=="ALIEN:1"],effects$prop.thr[effects$group=="ALIEN:0"])))
  
}

### frequency barplots with logarithmic scale
thresh.freq=function(effects=effects.glmSR.grass, data=species, y=T,ylim=c(0,100), leg = F) 
{
  tmp=effects
  
  # loop on the two groups : Native and Alien targets
  for (i in 1:2) {
    x=t(as.matrix(tmp[tmp$group==levels(tmp$group)[i],c("freq.thr", "freq.impact", "nb.sp","n.target" )]))
    x = log(x)
    barplot(x[3,]+1, names= paste("c",2:6, sep=""),las=1, ylim=c(0, max(log(ylim + 1))),
            col=c( "white"), border="grey", yaxt="n")
    barplot(x[4,]+1, col=c( "lightgrey"), border="grey",names="",axes=F,add=T)
    barplot(x[2,]+1, col=c( "grey"), names="",axes=F,add=T)
    barplot(x[1,]+1, col=c("black"),  names="",axes=F, add=T)
    axis(side=2, at=log(c(1,2,5,10,25,50,100, 200))+ 1, label=c(1,2,5,10,25,50,100,200), las=1)
    box(bty="l")
    if (y==T) mtext(side=2, text="Number of species", line=3.4, cex=0.8)
    
    if (leg == T & i==1) {
      legend('topright', bty="0", bg="white",legend=c("threshold", "significant", "ns. occurrences", "all species"),
             fill=c("black", "grey", "lightgrey", "white"), border= c("black", "black", "grey", "grey"), cex=0.9)
    }
  }
  
  chisq.test(t(cbind(effects$freq.thr[effects$group=="ALIEN:1"],effects$freq.thr[effects$group=="ALIEN:0"])))      
}


#####  Plotting impact size and prevalence
plot.impact = function(x="prevalence", y="nb.plot.impact", square = F){
  par(mfrow=c(1,3), cex=0.8, oma=c(2,2,2,2))
  
  for (i in 1:3) {
    M <- list(glmSR.grass, glmSRnat.grass, glmSRali.grass) [[i]]
    sp <- which(rownames(M$thresh)%in%aliens)
    M=lapply(M, FUN=function(x) x=x[sp,])
    
    if (square==T)  ylim <- xlim <- c(0, max(c(M$thresh[,x],M$thresh[,y]), na.rm=T))
    if (square==F)  ylim <- c(0, max(M$thresh[,y], na.rm=T)) ;  xlim <- c(0, max(M$thresh[,x], na.rm=T))
    
    plot(M$thresh[,x],M$thresh[,y] , ann=F,pch=20,  cex=1.2, xlim = xlim, ylim= ylim,
         col=c("forestgreen", "firebrick")[ rownames(M$thresh)%in% aliens +1] ) 
    if (square==T) abline(0,1)
    
    title( main= c("Total richness", "Native Richness","Alien Richness")[i])
  }
  mtext(1, text = x, outer= T, line=1)
  mtext(2, text = y, outer= T, line=1)
}


#####  Plotting impact size and prevalence
rank.impact = function(x="wtd.mean.diff", y="nb.plot.impact"){
  par(mfrow=c(1,3), cex=0.8, oma=c(2,2,2,2))
  for (i in 1:3) {
    M <- list(glmSR.grass, glmSRnat.grass, glmSRali.grass) [[i]]
    sp <- which(rownames(M$thresh)%in%aliens & ( M$thresh$nb.plot.impact>0))
    M=lapply(M, FUN=function(x) x=x[sp,])
    
    M$thresh$impact <- M$thresh[,x] * M$thresh[,y]
    M$thresh <- M$thresh[order(M$thresh$impact, decreasing=T),]
    
    plot(1:length(sp), M$thresh$impact,ann=F,type="h",cex=1.2,
         ylim= c(0,max(M$thresh$impact)),
         col=c("forestgreen", "firebrick")[ rownames(M$thresh)%in% aliens +1], xaxt="n") 
    axis(side=1, at=1:length(sp), label=rownames(M$thresh), las=3, cex.axis=0.7)
    title( main= c("Total richness", "Native Richness","Alien Richness")[i])
  }
  mtext(2, text = paste("impact index =", x, "*", y) , outer= T, line=1)
  mtext(2, text = paste("impact index =", x, "*", y) , outer= T, line=1)
}
=======
### Threshold analysis plotting functions

## plot summary outputs
plot.effect.summary <- function(  effects = glmSR.sum$class.summary ) {       
  
  n <- length(unique(effects$group)) 
  
  par(mfcol=c(n,3), oma=c(2,10,4,0), mar=c(2,2,1,1))
  for (j in 1:3) {
    effects = list(glmSR.sum$class.summary,glmSRnat.sum$class.summary ,glmSRali.sum$class.summary )[[j]] 
  for (i in 1:n) {
    S <- as.data.frame(effects[effects$group == unique(effects$group)[i],])
    barplot(S$prop.impact*100, ylim=c(0,50), col="lightgrey")
    par(new=T)
    barplot(S$prop.thr*100, col="black" , ylim=c(0,50))
     
  }   
  }
#    mtext(text="Abundance class", side=1, outer=T, line=0.5)
   mtext(text="% observed species", side=2, outer=T, line=8)  
#   mtext(side=2, text=c("Native\ntargets","Alien\ntargets"), at =c (0.25, 0.75),
#         outer=T,adj=0.5,line=3, cex=0.8, las = 1)
  
  mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=4, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
  mtext(side=3, text=c("Total\nrichness","Native\nrichness", "Alien\nrichness"),line=0, at=c(0.18,0.52,0.85),adj=0.5, outer=T, las=1)
  mtext(text="Abundance class", side=1, outer=T, line=1)
  
legend('topright', bty="0", bg="white",legend=c("negative", "critical"),
fill=c("lightgrey", "black"), border= c("black", "black"), cex=0.9)

}

### plotting individual species
plot.sp.glm=function(sp="ACHMIL", M=glmSR, var="SR", db=db, bst =T, boxplots =T ) {  
  
  Z= db[db$SpeciesCode==sp,]
  
  if (bst ==T) M$thresh <-M$boot.thresh
  
  #alt boxplot
  (if(boxplots) {
  col=rep("white",6)
  col[ as.numeric( M$thresh[sp,"th"])]="red"  # identify the threshold
  boxplot(as.formula(paste(var, " ~ abun")), data=Z, xlim=c(0.5,6.5), outline=F,
          varwidth=T, col=col,  xaxt = "n", yaxt="n", ann=F)
  abline(h = median(Z[Z$abun == 1, var], na.rm=T), lty="dotted")
  }
  else {
     # alt points
    col = rep("white",length(Z$abun))
    col[ Z$abun == as.numeric( M$thresh[sp,"th"])]="black"  # identify the threshold
    Z$abun.jit = jitter(Z$abun, amount = 0.4)
    plot(as.formula(paste(var, " ~ abun.jit")), data=Z, xlim=c(0.5,6.5),
         type="n",xaxt = "n", yaxt="n", ann=F)
    abline(h = mean(Z[Z$abun == 1, var], na.rm=T), lty="dotted")
    points(as.formula(paste(var, " ~ abun.jit")), data=Z, xlim=c(0.5,6.5),
       col="black", bg=col,pch=21,  xaxt = "n", yaxt="n", ann=F, cex = 0.8)
  })

  axis(1,at = 1:6, label = names(M$n.obs), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  title(main=sp, cex.main=0.8, line=0.2)
 

  ## adding significantce per class
  if (bst ==T) Pclass <- M$boot[sp,grep("Pnegative", names(M$boot))]
  if (bst !=T) Pclass <- M$P[sp,]
#   mtext(side=3, line=-1.1, at=2:6,text=sapply( Pclass, FUN=p2star)[2:6], cex=0.8)
  
#   mtext(side=3, text=paste ("rho" ,round(M$spearman[sp,"rho"],2), p2star(M$spearman[sp,"p.val"])), adj=1, cex=0.7)
  
  return(M$thresh[sp,"th"])
}


### plotting Effect sizes for individual species 
plot.sp.ES=function(sp="ACHMIL", M=glmSR, var="SR", db=databp) {  
  
  Z <- db[db$SpeciesCode==sp,]
  es <- as.numeric(M$est[sp,])
  low <- as.numeric(M$boot[sp,grep("q2.5", names(M$boot))])[2:6]
  hi <- as.numeric(M$boot[sp,grep("q97.5", names(M$boot))])[2:6]
  
  n <- as.numeric(M$n.obs[sp,])[2:6]
   ## correct for small sample sizes
  small <- es
  small[which(n>=3)] <- NA
  es[which(n<3)] <- NA
  low[n<3] <- NA
  hi[n<3] <- NA
  
  lims <-  max(c(abs(low),abs(hi)), na.rm=T) +0.02
  
  
  col <- rep("white",6)
  col[ as.numeric( M$boot.thresh[sp,"th"])] <- "black"  # identify the threshold
  
  ## alternatively:
  plot(1:5,  rep(0,5),ylim = c(-lims,+lims),xlim = c(0.5, 5.5), type= "n", xaxt = "n", yaxt="n", ann=F)
  axis(1,at = 1:5, label = names(M$est), tcl= 0.1,mgp=c(1,0.5,0),las=1)
  axis(2, tcl= 0.1,  mgp=c(1,0.5,0), las=1)
  title(main=sp, cex.main=0.8, line=0.2)
  
  abline(h=0,lty="dotted")

  # if (!all(is.na(small))) points(1:5, small, pch=".")  

  ## Add bootstrapped CI
  arrows(1:5,low,   1:5,hi, lwd=1, code =3, length=0.05, angle=90)
  points(1:5,es, bg = col[2:6], pch=21)
  
  mtext(side=3, line=-1.2, at=1:5,text=n, cex=0.7, las=1, hadj=1)
  
    ## adding significantce per class
#   Pclass <- M$P[sp,]
#   mtext(side=3, line=-1.1, at=1:5,text=sapply( Pclass, FUN=p2star), cex=0.8)
#   mtext(side=3, text=paste ("rho" ,round(M$spearman[sp,"rho"],2), p2star(M$spearman[sp,"p.val"])), adj=1, cex=0.7)
  
  return(M$boot.thresh[sp,"th"])
}

## plotting all significant species
plot.glm=function(M=glmSR, var="SR", db=databp, sel.criteria = c("th.exist"), bst =T, ES = T,
                  panels = c(4,4), boxplots = T) {
  
  par(mfrow=panels, mar=c(2,2,2,1), oma=c(3,3,1,1))
  if (bst ==T)  M$thresh <- M$boot.thresh
  
  # select only species according to seletion criteria :
  if (sel.criteria == "spear") sel <- rownames(M$spearman) [M$spearman[,"p.val"] <=0.05 &  M$spearman[,"rho"] <0 ]
  if (sel.criteria == "dev1") sel <-  rownames(M$glms) [M$glms[,"dev.ratio"] > 0.05 &  !is.na(M$thresh[,"th"]) ]
  if (sel.criteria == "none") sel <- TRUE
  if (sel.criteria == "th.exist") sel <-  rownames(M$thresh) [ !is.na(M$thresh[,"th"]) ]
  
  ### Loop on selected species 
  for (sp in sel)  {
    
    if (ES == F) th=plot.sp.glm(sp= sp, M=M, var=var, db=db, bst = bst, boxplots =boxplots )
    if (ES == T) th=plot.sp.ES(sp= sp, M=M, var=var, db=db)
    
    print(paste(sp,":",th))
  }
}

### thredhold barplots
thresh.prop=function(effects=effects.glmSR, data=species, y=T,ylim=c(0,0.35)) {
  #graphical representation
  barplot(effects$prop.thr[effects$group=="ALIEN:0"], names.arg= paste("c",2:6, sep=""), 
          col="grey", ylim=ylim, las=1 )
  if (y==T) mtext(side=2, text="Prop. threshold effects",line=3.4, cex=0.8)
  barplot( effects$prop.thr[effects$group=="ALIEN:1"], names.arg= paste("c",2:6, sep=""),
           col="black" , ylim=ylim , las=1)
  if (y==T) mtext(side=2, text="Prop. threshold effects", line=3.4, cex=0.8)
  
  print(wilcox.test(effects$prop.thr[effects$group=="ALIEN:1"],effects$prop.thr[effects$group=="ALIEN:0"], paired=T))
  print(friedman.test(cbind(effects$prop.thr[effects$group=="ALIEN:1"],effects$prop.thr[effects$group=="ALIEN:0"])))
  
}

### frequency barplots with logarithmic scale
thresh.freq=function(effects=effects.glmSR, data=species, y=T,ylim=c(0,100), leg = F) 
{
  tmp=effects
  
  # loop on the two groups : Native and Alien targets
  for (i in 1:2) {
    x=t(as.matrix(tmp[tmp$group==levels(tmp$group)[i],c("freq.thr", "freq.impact", "nb.sp","n.target" )]))
    x = log(x)
    barplot(x[3,]+1, names= paste("c",2:6, sep=""),las=1, ylim=c(0, max(log(ylim + 1))),
            col=c( "white"), border="grey", yaxt="n")
    barplot(x[4,]+1, col=c( "lightgrey"), border="grey",names="",axes=F,add=T)
    barplot(x[2,]+1, col=c( "grey"), names="",axes=F,add=T)
    barplot(x[1,]+1, col=c("black"),  names="",axes=F, add=T)
    axis(side=2, at=log(c(1,2,5,10,25,50,100, 200))+ 1, label=c(1,2,5,10,25,50,100,200), las=1)
    box(bty="l")
    if (y==T) mtext(side=2, text="Number of species", line=3.4, cex=0.8)
    
    if (leg == T & i==1) {
      legend('topright', bty="0", bg="white",legend=c("threshold", "significant", "ns. occurrences", "all species"),
             fill=c("black", "grey", "lightgrey", "white"), border= c("black", "black", "grey", "grey"), cex=0.9)
    }
  }
  
  chisq.test(t(cbind(effects$freq.thr[effects$group=="ALIEN:1"],effects$freq.thr[effects$group=="ALIEN:0"])))      
}


#####  Plotting impact size and prevalence
plot.impact = function(x="prevalence", y=c("mean.dif"), square = F, prop ="", bst = T,xlab = x, ylab =y ){
  par(mfrow=c(1,3), cex=0.8, oma=c(2,2,2,2))
  
  for (i in 1:3) {
    M <- list(glmSR, glmSRnat, glmSRali) [[i]]
    var <- c("SR", "SRnat", "SRali")[i]
    #     sp <- which(rownames(M$thresh)%in%aliens)  ## selects only aliens
    #     M=lapply(M, FUN=function(x) x=x[sp,])
    if (bst) M$thresh <- M$boot.thresh

  M$thresh$Y <- (if (length(y)==3) {
  apply(M$thresh, 1, FUN=function(k){
    p <- parse(text = eval(paste( k[y[1]] , y[2] , k[y[3]]))) 
    eval(p)
  })
}
else {
if ("y" %in% prop){
  calc.prop.index(var="SR" , M = M, index = y)
}
else  M$thresh[,y]
}
)


M$thresh$X <- (if (length(x)==3) {
  apply(M$thresh, 1, FUN=function(k){
    p <- parse(text = eval(paste( k[y[1]] , y[2] , k[y[3]]))) 
    eval(p)
  })
}
else {
if ("x" %in% prop){
  calc.prop.index(var="SR" , M = M, index = x)
}
else M$thresh[,x]
}
)


if (square==T)  ylim <- xlim <- c(0, max(c(M$thresh[,"X"],M$thresh[,"Y"]), na.rm=T))
if (square==F)  ylim <- c(0, max(M$thresh[,"Y"], na.rm=T)) ;  xlim <- c(0, max(M$thresh[,"X"], na.rm=T))


sp <- which(M$thresh[,"Y"] == 0)
(if (length(sp)>0 ) {
    plot(M$thresh[-sp,"X"],M$thresh[-sp,"Y"] , ann=F,pch=21,  cex=1, xlim = xlim, ylim= ylim,
         col=c("forestgreen", "firebrick")[ rownames(M$thresh[-sp,])%in% aliens +1],
         bg=c("forestgreen", "firebrick")[ rownames(M$thresh[-sp,])%in% aliens +1]) 
    points(M$thresh[sp,"X"],M$thresh[sp,"Y"] , ann=F,pch=21,  cex=1,
           col=c("forestgreen", "firebrick")[ rownames(M$thresh[sp,])%in% aliens +1])
}
else {
  plot(M$thresh[,"X"],M$thresh[,"Y"] , ann=F,pch=21,  cex=1, xlim = xlim, ylim= ylim,
       col=c("forestgreen", "firebrick")[ rownames(M$thresh[,])%in% aliens +1],
       bg=c("forestgreen", "firebrick")[ rownames(M$thresh[,])%in% aliens +1]) 
})

    if (square==T) abline(0,1)
    
    title( main= c("Total richness", "Native Richness","Alien Richness")[i])
  }
  mtext(1, text = xlab, outer= T, line=0)
  mtext(2, text = ylab, outer= T, line=0)
}


#####  Plotting impact size and prevalence
rank.impact = function(x="index"){
  par(mfrow=c(1,3), cex=0.8, oma=c(2,2,2,2))
  for (i in 1:3) {
    M <- list(glmSR, glmSRnat, glmSRali) [[i]]
     M$thresh <- M$boot.thresh
    
    sp <- which(rownames(M$thresh)%in%aliens & ( M$thresh$nb.plot.impact>0))
    M=lapply(M, FUN=function(x) x=x[sp,])
    M$thresh$impact <-  M$thresh[,x]
    
    o  <- order(M$thresh$impact, decreasing=T)
    M=lapply(M, FUN=function(x) x=x[o,])
    
    t = seq(0,1, 0.01)
    S <- matrix(t, length(x),length(x), byrow=T)
    P <- matrix(t, length(x),length(x),byrow=F)
    z2 <-  S*P
    
    library(lattice)
    
    x11() 

    levelplot(z2~ S * P, col.regions= colorRampPalette(c("beige" , "firebrick", "purple")),
              main = "impact index" )
    trellis.focus("panel", 1, 1, highlight=FALSE)
    lpoints(M$thresh$prop.wt.mean.dif,M$thresh$prop.plot.impact , pch=3, col="black", cex=2)
    trellis.unfocus()
  
    x11()

    plot(1:length(sp), M$thresh$impact,ann=F,type="h",cex=1.2,
         ylim= c(0,max(M$thresh$impact)),
         col=c("forestgreen", "firebrick")[ rownames(M$thresh)%in% aliens +1], xaxt="n") 
    axis(side=1, at=1:length(sp), label=rownames(M$thresh), las=3, cex.axis=0.7)
    title( main= c("Total richness", "Native Richness","Alien Richness")[i])
  }
  mtext(2, text = paste("impact index =", x, "*", y) , outer= T, line=1)
  mtext(2, text = paste("impact index =", x, "*", y) , outer= T, line=1)
}
>>>>>>> FETCH_HEAD
