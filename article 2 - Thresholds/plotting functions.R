### Threshold analysis plotting functions

## plot summary outputs
plot.effect.summary <- function(  effects = effects.glmSR.grass ) {       
  
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
  
  ## adding significantce per class
  if (type == "boot") Pclass <- M$boot[sp,grep("Pnegative", names(M$boot))]
  if (type != "boot") Pclass <- M$P[sp,]
  mtext(side=3, line=-1.1, at=2:6,text=sapply(  Pclass, FUN=p2star)[2:6], cex=0.8)
  
  
 mtext(side=3, text=paste ("rho" ,round(M$spearman[sp,"rho"],2), p2star(M$spearman[sp,"p.val"])), adj=1, cex=0.7)
  mtext(side=3, text=paste ("rho" ,round(M$spearman[sp,"rho"],2), p2star(M$spearman[sp,"p.val"])), adj=1, cex=0.7)
  
  title(main=sp, cex.main=0.8, line=1)
  
  return(M$thresh[sp,"th"])
}

## plotting all significant species
plot.glm=function(M=glmSR, var="SR", db=databp, sel.criteria = c("spear"), type="boot") {
  
  par(mfrow=c(5,5), mar=c(2,2,2,2))
  if (type=='boot')  M$thresh <- M$boot.thresh
  
  # select only species according to seletion criteria :
  if (sel.criteria == "spear") sel <- rownames(M$spearman) [M$spearman[,"p.val"] <=0.05 &  M$spearman[,"rho"] <0 ]
  if (sel.criteria == "dev1") sel <-  rownames(M$glms) [M$glms[,"dev.ratio"] > 0.05 &  !is.na(M$thresh[,"th"]) ]
  if (sel.criteria == "none") sel <- TRUE
  if (sel.criteria == "th.exist") sel <-  rownames(M$thresh) [ !is.na(M$thresh[,"th"]) ]
  
  ### Loop on selected species 
  for (sp in sel)  {
    th=plot.sp.glm(sp= sp, M=M, var=var, db=db, type = type)
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
#     sp <- which(rownames(M$thresh)%in%aliens)  ## selects only aliens
#     M=lapply(M, FUN=function(x) x=x[sp,])
    
    if (square==T)  ylim <- xlim <- c(0, max(c(M$thresh[,x],M$thresh[,y]), na.rm=T))
    if (square==F)  ylim <- c(0, max(M$thresh[,y], na.rm=T)) ;  xlim <- c(0, max(M$thresh[,x], na.rm=T))
    
sp <- which(M$thresh[,y] == 0)
(if (length(sp)>0 ) {
    plot(M$thresh[-sp,x],M$thresh[-sp,y] , ann=F,pch=21,  cex=1, xlim = xlim, ylim= ylim,
         col=c("forestgreen", "firebrick")[ rownames(M$thresh[-sp,])%in% aliens +1],
         bg=c("forestgreen", "firebrick")[ rownames(M$thresh[-sp,])%in% aliens +1]) 
    points(M$thresh[sp,x],M$thresh[sp,y] , ann=F,pch=21,  cex=1,
           col=c("forestgreen", "firebrick")[ rownames(M$thresh[sp,])%in% aliens +1])
}
else {
  plot(M$thresh[,x],M$thresh[,y] , ann=F,pch=21,  cex=1, xlim = xlim, ylim= ylim,
       col=c("forestgreen", "firebrick")[ rownames(M$thresh[,])%in% aliens +1],
       bg=c("forestgreen", "firebrick")[ rownames(M$thresh[,])%in% aliens +1]) 
})

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
