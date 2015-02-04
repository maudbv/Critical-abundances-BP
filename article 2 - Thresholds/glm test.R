### GLM on SR ~ abund analyses

#### functions

# fitting GLMs
#  = testing difference in a community metric (e.g. SR) between low abundances
# and increasing abundance classes

glm.test=function(db=databp, var='SR',domclass='abun',zeros=F, verbose=F, test="w",
                  env=envplot, min.occur=10,  min.class=3, alpha=0.05) {
  require(coin)
    
   ## sp spanning the minimum number of abundance classes
   a <- unique(names(which(rowSums(table(db$SpeciesCode, db$abun)>0)>=min.class)))  
 
  ## species with minimum occurrences
    # calculate species prevalence
    db$Sp.occurence <- table(as.character(db$SpeciesCode))[match(as.character(db$SpeciesCode),unique(as.character(db$SpeciesCode)))]
   b <- unique(names(which(db$Sp.occurence>min.occur)))
  
  # select only species in a and b groups :
  db.modif <- db[which( (db$SpeciesCode %in% a) & (db$SpeciesCode %in% b)),] 
  
  # list of species name to be included in analysis :
  sp.names <- unique(db.modif$SpeciesCode)
  
  # initiate results
  init <-  data.frame(matrix(NA, nrow=length(sp.names), ncol= 5,
                           dimnames=list(sp.names,c("c2", "c3", "c4", "c5", "c6"))))
  diff <-  n.obs <- est<- P <- z <- init
  
  glms <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 5,
                           dimnames=list(sp.names, c("min.sig" ,"th", "df", "resid.dev","dev.ratio"))))
  
  spear <-data.frame(matrix(NA, nrow=length(sp.names), ncol= 2,
                            dimnames=list(sp.names, c("rho" ,"p.val"))))
                     
  
  # looping on species
  for (i in 1:length(sp.names) ) {
    
    sp=sp.names[i]  # select species name
    x=db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",var,'PlotName')] # select occurrences of the species
    names(x)=c('abun','var','PlotName')
    
    #    # add all other plots as absences :
    #     z= cbind(abun=NA,env[which(!env$PLOTID %in% x$PlotName),c(var,"PLOTID")]) 
    #     x=apply(x,2,as.numeric)
    #     z=apply(z,2,as.numeric)
    #     x=as.data.frame(rbind(x,z))
    
abun=sort(as.numeric(as.character(na.omit(unique(x$abun))))) ## list of abundance classes for species i
n.obs[i,]=sapply(2:6, FUN=function(j) length(x[which(x$abun==j),"var"]))
    
         
      #calculate difference in class mean SR
      diff[i,] =  sapply(2:6, FUN=function(j) mean(x[which(x$abun==j),"var"], na.rm=T) - mean (x[which(x$abun==min(abun)),"var"], na.rm=T))
      
      # spearmnan test across all classes
      s =cor.test(x$var,x$abun, method="spearman")  
      spear[i,] =  c(s$estimate, s$p.value)
      
      # GLM test
      f <-  glm(x$var ~ as.factor(x$abun), family=poisson(log))
      
      # detect thresholds
      th=NA
      min.sig=NA
      sig  <-  which((summary(f)$coef [, 4] <= alpha) & ( summary(f)$coef [, 3] <0)) ## which classes show signif negative diff
      neg <-   which(summary(f)$coef [, 3]<0)  # which classes show negative diff
      
      if (length(sig)!=0) {
        min.sig = min(sig)    # lowest class with signif negative difference
        test=F
        j=0
        
        while (test==F & j<=5) {   # search for thresholds in case min.sig is not followed by negative differences
          if (all(((min.sig+j):max(abun)) %in% neg)) { th=min.sig ; test=T }
          j=j+1                  
        }
      }
      
      glms[i,] <-  c(min.sig = min.sig, th= th, df= f$df.resid, resid.dev= f$dev,
                     dev.ratio= (f$null.deviance -f$dev)/f$null.deviance  )
      
      n = 1:(length(abun)-1)
      est[i,n] <- summary(f)$coef [-1, 1]
      z[i,n] <- summary(f)$coef [-1, 3]
      P[i,n] <- summary(f)$coef [-1, 4]
      
    }      
   
return(list(glms=glms, spearman=spear, n.obs = n.obs,  diff = diff,est= est, z=z,P= P))
}


# summarizing results
summary.glmtest = function(M = glmSR,data=species, group="ALIEN", graph=T) {
  data=data[rownames(M$diff),]
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
  
  par(mfrow=c(G,2), oma=c(4,6,6,2), mar=c(2,2,1,1))
  out=NULL   
  
  for (j in 1:G) {
    S=sub[[j]]
    
    # number of potentieally signif species occurrence per class
    n.sp=sapply(names(S$diff), FUN= function(x) {
      length(which(!is.na(S$diff[,x])))
    })
    
    # number of significant negative effects per class
    n.impact=sapply(names(S$diff), FUN= function(x) {
      length(which(S$P[,x]<=0.05 & S$z[,x]<0))
    })
    
    # number of significant positive effects per class
    n.positive=sapply(names(S$diff), FUN= function(x) {
      length(which(S$P[,x]<0.05 & S$z[,x]>=0))
    })
    
    
    # proportion of signifi negative effects per class
    p.impact=n.impact / n.sp
    
    # number of times each class is the threshold
    
    freq.thr =table(as.factor(S$glms$th))[match(2:6,names(table(as.factor(S$glms$th))))]
    names(freq.thr)=names(n.sp)
    freq.thr[is.na(freq.thr)]=0
    prop.thr = freq.thr / n.sp
    if(sum(freq.thr, na.rm=T)!=0) perc.thr = freq.thr / sum(freq.thr, na.rm=T)
    if(sum(freq.thr, na.rm=T)==0) perc.thr = freq.thr
    
    #graphical representation
    if (graph) {
      
      barplot(p.impact)
      if (j==1) mtext(side=3, text="% negative\neffects",adj=0.5,line=1, cex=0.8)
      
      barplot(prop.thr, col="black" )
      if (j==1) mtext(side=3, text="% threshold\neffects", adj=0.5,line=1, cex=0.8)
      
      out=rbind(out, data.frame(group=names [j],nb.sp=n.sp, freq.impact=n.impact, freq.positives =n.positives,
                                prop.impact=p.impact,freq.thr=freq.thr, prop.thr=prop.thr, perc.thr=perc.thr))
    }
    
    
    mtext(side=2, text=c("Alien\ntargets", "Native\ntargets"),line=3, at=c(0.3,0.8),adj=0.5, outer=T, las=1)
    mtext(text="Abundance class", side=1, outer=T, line=2)
  }
  return(out)
}


### plotting individual species
plot.sp.glm=function(sp="ACHMIL", M=glmSR, var="SR", db=databp)
{  
  Z= db[db$SpeciesCode==sp,]
  col=rep("white",6)
  col[ as.numeric( M$glms[sp,"th"])]="red"
  boxplot(as.formula(paste(var, " ~ abun")), data=Z, xlim=c(0,6), outline=F, varwidth=T, col=col)
  mtext(side=3, line=-1.1, at=2:6,text=sapply(M$P[sp,], FUN=p2star), cex=0.8)
  mtext(side=3, text=paste ("rho" ,round(M$spearman[sp,"rho"],2), p2star(M$spearman[sp,"p.val"])), adj=1, cex=0.7)
  title(main=sp, cex.main=0.8)

  return(M$glms[sp,"th"])
}

## plotting results
plot.glm=function(M=glmSR, var="SR", db=databp, sel.criteria = c("spear")) {
  par(mfrow=c(4,5), mar=c(2,2,2,2))
  for (sp in rownames(M$diff)){
    # only species with significant kruskal walis test :
    if (sel.criteria == "spear") sel <- M$spearman[sp,"p.val"] <=0.05 &  M$spearman[sp,"rho"] <0
#     if (sel.criteria == "dev1") sel <- M$glms[sp"dev.ratio"] > 0.05 &  !is.na(M$glms[sp,"th"])
    if (sel.criteria == "none") sel <- TRUE
    
    if (sel==T){
      th=plot.sp.glm(sp= sp, M=M, var=var, db=db)
      print(paste(sp,":",th))
    }
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


### logarithmic scale

thresh.freq=function(effects=effects.glmSR, data=species, y=T,ylim=c(0,100), leg=F) {
  tmp=effects
  # Native targets
  x=t(as.matrix(tmp[tmp$group=="ALIEN:0",c("freq.thr", "freq.impact", "nb.sp")]))
  #    x[x!=0] = log(x[x!=0])
  x = log(x)
  barplot(x[3,]+1, names= paste("c",2:6, sep=""),las=1, ylim=c(0, max(log(ylim + 1))),
          col=c( "white"), yaxt="n")
  barplot(x[2,]+1, col=c( "grey"), names="",axes=F,add=T)
  barplot(x[1,]+1, col=c("black"),  names="",axes=F, add=T)
  axis(side=2, at=log(c(1,2,5,10,25,50,100, 200))+ 1, label=c(1,2,5,10,25,50,100,200), las=1)
  box(bty="l")
  if (y==T) mtext(side=2, text="Number of species", line=3.4, cex=0.8)
  
  if (leg) {
    legend('topright', bty="0", bg="white",legend=c("threshold", "impact", "no impact"),
           fill=c("black", "grey", "white"), cex=0.9)
  }
  
  ### Alien targets
  x=t(as.matrix(tmp[tmp$group=="ALIEN:1",c("freq.thr", "freq.impact", "nb.sp")]))
  #    x[x!=0] = log(x[x!=0])
  x = log(x)
  barplot(x[3,]+1, names= paste("c",2:6, sep=""),las=1, ylim=c(0, max(log(ylim + 1))),
          col=c( "white"), yaxt="n")
  barplot(x[2,]+1, col=c( "grey"), names="",axes=F,add=T)
  barplot(x[1,]+1, col=c("black"),  names="",axes=F, add=T)
  axis(side=2, at=log(c(1,2,5,10,25,50,100, 200))+ 1, label=c(1,2,5,10,25,50,100,200), las=1)
  box(bty="l")
  if (y==T) mtext(side=2, text="Number of species", line=3.4, cex=0.8)
  
  chisq.test(t(cbind(effects$freq.thr[effects$group=="ALIEN:1"],effects$freq.thr[effects$group=="ALIEN:0"])))      
}

