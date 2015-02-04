# null model: assign a random abundance to species within each comm

simul1.class.rhos=function(db=databp,var="SR",status.shuffle=T, status='ALIEN', out='rhos', recalculate=T,min=10,alpha=0.01, nreps=99)  {

require(doBy)

  db$Sp.occurence=table(as.character(db$SpeciesCode))[match(as.character(db$SpeciesCode),
                                                            unique(as.character(db$SpeciesCode)))]
  db=db[db$Sp.occurence>min,]

  output=data.frame(matrix(NA, ncol=nreps+1, nrow=length(unique(db$SpeciesCode))))

# OBSERVED values
    classSR=classSR.fun(db=db, var=var, min.occur=min, min.class=3, alpha=alpha)
    output[,1]= classSR$rho
  
## SIMULATIONS  
  for (i in 1:nreps) {
    
    #Null shuffling of abundance classes within sites
    nul=tapply(as.character(db$domlevels),as.character(db$PlotName), FUN=function(x) sample(x))
    nul.db=db
    nul.db$domlevels=as.numeric(as.character(unlist(nul)))
  
    if (recalculate==T) {
    als= sapply(nul,FUN=function(x) sum(as.numeric(x), na.rm=T))
    names(als)=unique(as.character(db$PlotName))
    nul.db$SRali= als[as.character(nul.db$PlotName)]
    
    nats= sapply(nul,FUN=function(x) sum(as.numeric(x), na.rm=T))
    names(nats)=unique(as.character(db$PlotName))
    nul.db$SRnat= nats[as.character(nul.db$PlotName)]
    }
    
    tmp=nul.db
    null.classSR=classSR.fun(db=tmp, var=var,  min.occur=min, min.class=3, alpha=0.01)
    output[,i+1]=  null.classSR$rho
}   

  rownames(output)= unique(db$SpeciesCode) 
  
  return(output)
}


results.simul1.rhos=function(x=sim1.rhos) {
  x=data.frame(rho= -x[,1], mean.null= - rowMeans(x, na.rm=T))
  x$ALIEN= species[rownames(x), "ALIEN"]
  
  plot(c(0.5,2.5),c(min(x[,1], na.rm=T),max(x[,1], na.rm=T)), type='n', ann=F, axes=F,ylim=c(-0.8,0.8))
  abline(h=0)
  
  quant=cbind(quantile(x[x$ALIEN==1,1], probs=c(0.0275,0.9725), na.rm=T),quantile(x[x$ALIEN==0,1], probs=c(0.0275,0.9725), na.rm=T))
  arrows(1:2, quant[1,],1:2, quant[2,], col='black',angle = 90,code =3, length=0.05)  

  tmp=boxplot(rho ~ ALIEN , x, varwidth=T, ann=F,axes=F,col=c('lightgrey', 'grey50'),
           pars = list(boxwex = 0.5, staplewex = 0, outwex = 0),
         whisklty='blank',  staplelty='blank', boxlty='solid',outline=F,add=T,
          names=NULL,xaxt='n', cex.axis=0.8, las=1)
  
  out=rbind(x[which(x$ALIEN==0 & (x$rho<quant[1,1] | x$rho>quant[2,1])),],
            x[which(x$ALIEN==1 & (x$rho<quant[1,2] | x$rho>quant[2,2])),])

  points(out$ALIEN+1,out$rho, pch='.', cex=2)
  axis(1,at=1:2, label=c('Natives', 'Aliens'), cex.axis=0.8, tcl=-0.3, mgp=c(3,0.3,0))
# boxplot(mean.null ~ ALIEN , x, varwidth=T, ann=F, axes=F,
#             col=NA,border='blue', ylim=c(-0.8,0.8), add=T)
box(bty='o')
  # difference from null model mean
  y= quantile(x$rho[x$ALIEN==0], 0.75, na.rm=T) +0.1
t=t.test(x$rho[x$ALIEN==0],x[x$ALIEN==0,2],alternative = "greater", paired=F)
  if(t$p.val>0.05) {
  t=t.test(x$rho[x$ALIEN==0],x[x$ALIEN==0,2],alternative = "less", paired=F)
  }
text(x=1.2, y=y, label=(p2star(t$p.val)), cex=0.7)
 

y= quantile(x$rho[x$ALIEN==1], 0.75, na.rm=T) +0.1
t=t.test(x$rho[x$ALIEN==1],x[x$ALIEN==1,2],alternative = "greater", paired=F)
  if(t$p.val>0.05) {
  t=t.test(x$rho[x$ALIEN==1],x[x$ALIEN==1,2],alternative = "less", paired=F)
  }
text(x=2.2, y=y, label=(p2star(t$p.val)), cex=0.7)
  
# # difference between alien and natives
# t=t.test(x$rho[x$ALIEN==1],x[x$ALIEN==0,1],alternative = "greater")
#   if(t$p.val<=0.05) p=(p2star(t$p.val))
# 
#   if(t$p.val>0.05) {
#   t=t.test(x$rho[x$ALIEN==1],x[x$ALIEN==0,1],alternative = "less")
#   p=(p2star(t$p.val))
#   }
#   
#   mtext(side=1,text=substitute(italic("t.test: ")*p,list(p=p)),line=-1, adj=1, cex=0.6)
}


