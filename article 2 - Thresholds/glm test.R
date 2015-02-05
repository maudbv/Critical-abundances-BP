### FUNCTIONS for GLM analyses of SR ~ abund classes

### ANALYSES 

# fitting GLMs
#  = testing difference in a community metric (e.g. SR) between low abundances
# and increasing abundance classes using Generalized linear models 

glm.test <- function(db=databp[databp$PlotName %in% realgrasslands,], var='SR',
                     min.occur = 4,  min.class = 3, alpha=0.05) {
  
  ## sp spanning the minimum number of abundance classes in which they have a minimum occurrence
  a <- unique(names(which(rowSums(table(db$SpeciesCode, db$abun)>=min.occur)>=min.class)))  
  
  # select only species in a and b groups :
  db.modif <- db[which(db$SpeciesCode %in% a),] 
  
  # list of species name to be included in analysis :
  sp.names <- unique(db.modif$SpeciesCode)
  
  # initiate results
  init <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 5,
                             dimnames=list(sp.names,c("c2", "c3", "c4", "c5", "c6"))))
  diff <-  n.obs <- est<- P <- z <- init
  
  glms <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 3,
                           dimnames=list(sp.names, c( "df", "resid.dev","dev.ratio"))))
  
  spear <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 2,
                            dimnames=list(sp.names, c("rho" ,"p.val"))))
  
  thresh <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 9,
                             dimnames=list(sp.names, c("min.sig" ,"th",
                                                       "prevalence", "nb.plot.impact","prop.plot.impact",
                                                       "mean.diff","wtd.mean.diff", "th.diff", "max.diff"))))
  
  # looping on species
  for (i in 1:length(sp.names) ) {
    
    sp <- sp.names[i]  # select species name
    x <- db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",var,'PlotName')] # select occurrences of the species
    names(x) <- c('abun','var','PlotName')
    
    abun <- sort(as.numeric(as.character(na.omit(unique(x$abun))))) ## list of abundance classes for species i
    n.obs[i,] <- sapply(2:6, FUN=function(j) length(x[which(x$abun==j),"var"]))
    
    
    #calculate difference in class mean SR
    diff[i,] <-  sapply(2:6, FUN=function(j) mean(x[which(x$abun==j),"var"], na.rm=T) - mean (x[which(x$abun==min(abun)),"var"], na.rm=T))
    
    # spearmnan test across all classes
    s <- cor.test(x$var,x$abun, method="spearman")  
    spear[i,] <-  c(s$estimate, s$p.value)
    
    # GLM test
    f <-  glm(x$var ~ as.factor(x$abun), family=poisson(log))
    glms[i,] <-  c(df= f$df.resid, resid.dev= f$dev,dev.ratio= (f$null.deviance -f$dev)/f$null.deviance  )
    
    n <-  1:(length(abun)-1)
    est[i,n] <- summary(f)$coef [-1, 1]
    z[i,n] <- summary(f)$coef [-1, 3]
    P[i,n] <- summary(f)$coef [-1, 4]
    
    #### Detect thresholds
    
    # which classes show negative diff
    neg <-   which(summary(f)$coef [, 3]<0) 
    
    ## which classes show signif negative diff
    sig  <-  which((summary(f)$coef [, 4] <= alpha) & ( summary(f)$coef [, 3] <0)) 
    
    #### Threshold detection
    thresh[i,] <- (function() {
    
    min.sig <- th <- prevalence <- n.plot.impact <- mean.diff <- wtd.mean.diff <- th.diff <- max.diff <- NA
    
    ### identifying a threshold 
    ###  => lowest abundance where significant negative impact, followed by only negative trends (significant or not).
    if (length(sig)!=0) {
      min.sig <-  min(sig)    # init threshold value with lowest class with signif negative difference
      test <-  F
      j <- 0
      
      # search for thresholds in case min.sig is not followed by negative differences
      th.sig = sapply(sig, FUN= function(k) {
        c1 <- all(((k):max(abun)) %in% neg) # all higher classes have negative differences
        c2 <- (length(x[x$abun == k, "var"]) >= min.occur) # at least min.occur observation in the threshold class
        return(c1 & c2)
      })
      
      if (! all(th.sig == F)) th = min(sig[th.sig == T])
    }     
    
    
    ## additional stats per species :
    prevalence <- dim(x)[1]
    n.plot.impact <- sum(x$abun >= th, na.rm=T)
    prop.plot.impact <- n.plot.impact/prevalence
    
    #impact size for species with a threshold of impact :
    if (!is.na(th))  { 
    
        # mean impact
        mean.diff <- -mean(as.numeric(diff[i,c(th:6)-1]), na.rm=T)
        
        # frequency weighted mean impact
        d <- as.numeric(diff[i, c(th:6)-1])
        nb <- as.numeric(n.obs[i, c(th:6)-1])
        wtd.mean.diff <- -sum(d*nb, na.rm=T)/sum(nb, na.rm=T)
        
        # Max impact
        #! maximum diff is in fact the minimum because negative values
        max.diff <- - min(as.numeric(diff[i, c(th:6)-1]), na.rm=T) 
  
        # Threshold difference :
        th.diff = -as.numeric(diff[i,th-1])
    }
    
  return( c(min.sig, th, prevalence, n.plot.impact, prop.plot.impact, mean.diff,  wtd.mean.diff,th.diff,max.diff))
    })()      
  }
  return(list(glms=glms, thresh = thresh, spearman=spear, n.obs = n.obs,  diff = diff,est= est, z=z,P= P))
}

# summarizing results on thresholds per group of species
summary.glmtest = function(M = glmSR.grass,data=species, group="ALIEN") {
  
  ### select onlys species in the model
  data=data[rownames(M$diff),]
  
  ## divide species according to grouping factor (ie Alien vs Native)
  if (is.null(group)) {
    sub=M
    G=1
  }
  
  if (!is.null(group)) { 
    sub=list()
    i=1
    G=length(unique(data[,group]))
    for (g in sort(unique(data[,group]))) {
      sub[[i]]=lapply(M, FUN= function(x) x[rownames(x) %in% rownames(data)[data[,group]==g],] )
      i=i+1
    }
    names=paste(group, ":",sort(unique(data[,group])), sep="")
  }                                 
  
  out=NULL   
  
  for (j in 1:G) {
    S=sub[[j]]
    
    # number of potentieally signif species occurrence per class
    n.sp=sapply(names(S$diff), FUN= function(x) {
      length(which(!is.na(S$diff[,x])))
    })
    
    # number of species which have a threshold per class
    n.target=sapply(names(S$diff), FUN= function(x) {
      length(which(!is.na(S$diff[,x]) & !is.na(S$thresh$th)) )
    })
    
    
    # number of significant negative effects per class
    n.negative=sapply(names(S$diff), FUN= function(x) {
      length(which(S$P[,x]<=0.05 & S$z[,x]<0))
    })
    
    
    # number of significant positive effects per class
    n.positive=sapply(names(S$diff), FUN= function(x) {
      length(which(S$P[,x]<0.05 & S$z[,x]>=0))
    })
    
    # number of significant negative impavts per class for species with a threshld
    n.impact=sapply(names(S$diff), FUN= function(x) {
      
      length(which(S$P[,x]<=0.05 & S$z[,x]<0 & !is.na(S$thresh$th)))
    })
    
    # proportion of signifi negative effects per class
    p.impact=n.impact/ n.sp
    
    # proportion of signifi negative effects per class
    p.negative=n.negative/ n.sp
    
    # number of times each class is the threshold
    freq.thr =table(as.factor(S$thresh$th))[match(2:6,names(table(as.factor(S$thresh$th))))]
    names(freq.thr)=names(n.sp)
    freq.thr[is.na(freq.thr)]=0
    prop.thr = freq.thr / n.sp
    if(sum(freq.thr, na.rm=T)!=0) perc.thr = freq.thr / sum(freq.thr, na.rm=T)
    if(sum(freq.thr, na.rm=T)==0) perc.thr = freq.thr 
    
    out=rbind(out, data.frame(group=names [j],nb.sp=n.sp,n.target = n.target, freq.impact=n.impact, freq.negative =n.negative,freq.positive =n.positive,
                              prop.impact=p.impact,freq.thr=freq.thr, prop.thr=prop.thr, perc.thr=perc.thr))
  }
  return(out)
}

### GRAPHS

## plot summary outputs
plot.effect.summary <- function(  effects = effects.glmSR ) {       
  
  n <- length(unique(effects$group)) 
  par(mfrow=c(n,2), oma=c(4,5,5,2), mar=c(2,2,1,1))
  
  for (i in 1:n) {
    S <- effects[effects$group == unique(effects$group)[i],]
    barplot(S$prop.impact*100, ylim=c(0,50))
    if (i==1) mtext(side=3, text="% negative\neffects",adj=0.5,line=1, cex=0.8)
    
    barplot(S$prop.thr*100, col="black" , ylim=c(0,50))
    if (i==1) mtext(side=3, text="% threshold\neffects", adj=0.5,line=1, cex=0.8)
    
  }   
  mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=3, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
  mtext(text="Abundance class", side=1, outer=T, line=2)
}

### plotting individual species
plot.sp.glm=function(sp="ACHMIL", M=glmSR, var="SR", db=databp) {  
  Z= db[db$SpeciesCode==sp,]
  col=rep("white",6)
  col[ as.numeric( M$thresh[sp,"th"])]="red"
  boxplot(as.formula(paste(var, " ~ abun")), data=Z, xlim=c(0,6), outline=F, varwidth=T, col=col)
  mtext(side=3, line=-1.1, at=2:6,text=sapply(M$P[sp,], FUN=p2star), cex=0.8)
  mtext(side=3, text=paste ("rho" ,round(M$spearman[sp,"rho"],2), p2star(M$spearman[sp,"p.val"])), adj=1, cex=0.7)
  title(main=sp, cex.main=0.8)
  
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
      legend('topright', bty="0", bg="white",legend=c("threshold", "impact", "no impact", "all species"),
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

