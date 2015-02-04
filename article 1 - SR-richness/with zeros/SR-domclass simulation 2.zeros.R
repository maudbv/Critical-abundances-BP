
################################################################################################################  

######## SIMULATION 2 : random distribtuion of alien/native status
# # 
shuff.all=shuff.classSR(db = databp,species = species, envplot=envplot,zeros=T,min.occur = 10, min.class = 3, alpha = 0.01, nreps = 999) ; save(shuff.all, file='shuff.all.Rdata')
shuff.grass=shuff.classSR(db = databp[databp$PlotName%in% grasslands,],envplot=envplot,zeros=T,species = species, min.occur = 10, min.class = 3, alpha = 0.01, nreps = 999) ; save(shuff.grass, file='shuff.grass.Rdata')
shuff.wood=shuff.classSR(db = databp[databp$PlotName%in% woodlands,],envplot=envplot,zeros=T,species = species, min.occur = 10, min.class = 3, alpha = 0.01, nreps = 999) ; save(shuff.wood, file='shuff.wood.Rdata')

# x11()
# results.shuff(shuff.table=shuff.all, ses='ses')
# results.shuff(shuff.table=shuff.grass, ses='ses')
# results.shuff(shuff.table=shuff.wood, ses='ses')

#  # By subcomponent
# par(mfrow=c(1,3), oma=c(1,3,1,1), mar=c(3,1,1,1), cex=0.9, las=1)
# x=sort(c(seq(1.2,12.2, 2), seq(1.8,12.8,2)))
# 
# names=c("Natives.negSR","Aliens.negSR","Natives.posSR","Aliens.posSR")
# y=data.frame(SR=shuff.all[1:1000,names], grass= shuff.grass[1:1000,names],wood=shuff.wood[1:1000,names])
#     # significant differences ?
#    sig.inf=sapply(1:12,FUN=function(x) {
#       p=(sum((y[,x]<y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
#       return(p)   
#      } )      
#                                          
#   sig.sup=sapply(1:12,FUN=function(x) {
#        p=(sum((y[,x]>y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
#       return( p )   
#      } )      
# # SES
# ses=sapply(1:12,FUN=function(x) {
#       ses= (y[1,x] - mean(y[,x], na.rm=T)) /sd(y[,x], na.rm=T)
#       return(ses)   
#      } )  
# 
# plot(x, ses, type='h', lwd=10, lend=1, col=c('grey', 'black'), ann=F, axes=F, xlim=c(1,12),
#      ylim=c(-4,13))
# abline(h=0)
# axis(1, at=c(1.5,3.5,5.5,7.5,9.5,11.5), line=-0.5, labels = rep(c('-','+'),3), tcl=-0.3, lty=0, cex.axis=1.2)
# axis(1, at=c(2.5,6.5,10.5), line=0.5,lty=0,hadj=0.5,labels = c('All plots','Grasslands', 'Woodlands'),cex.axis=0.8)
# axis(1, at=c(0.5,4.5,8.5,12.5), lwd=0,lwd.ticks=1, labels = NA)
# axis(1, at=c(2.5,6.5,10.5), lwd=0,lwd.ticks=1, labels = NA, tcl=-0.2)
# axis(2, las=1, cex.axis=0.8, mgp=c(3,0.7,0), tcl=-0.4 )
# # title(ylab='Standardized frequency', line=2, cex.lab=0.9)
# abline(v=4.5, lty='dotted')                                                    
# abline(v=8.5, lty='dotted')
# box(bty='o')  
# legend(x=8.5, y=10, legend=c('Natives','Aliens'), fill=c('grey', 'black'), cex=0.7, bty='n')
# 
# sig=apply(cbind(sig.inf, sig.sup),1, min, na.rm=T)
# sig=sapply(sig, p2star)
# sig[sig=='ns']=''
# text(y=ses[ses>0]+0.5, x=x[ses>0], label=sig[ses>0], cex.lab=0.7)
# text(y=ses[ses<0]-0.5, x=x[ses<0], label=sig[ses<0], cex.lab=0.7)
# mtext(side=1,line=-1, at=1:12, text=y[1,], cex=0.7, col='grey40', font=3)
# mtext(side=3,line=0,adj=0.5,  text="a) Total species richness", cex=0.9)
# 
# # Native richness
# names=c("Natives.negSRnat", "Aliens.negSRnat", "Natives.posSRnat" , "Aliens.posSRnat")
# y=data.frame(SR=shuff.all[1:1000,names], grass=shuff.grass[1:1000,names], wood=shuff.wood[1:1000,names])
#     # significant differences ?
#    sig.inf=sapply(1:12,FUN=function(x) {
#       p=(sum((y[,x]<y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
#       return(p)   
#      } )      
#                                          
#   sig.sup=sapply(1:12,FUN=function(x) {
#        p=(sum((y[,x]>y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
#       return( p )   
#      } )      
#  
# x=sort(c(seq(1.2,12.2, 2), seq(1.8,12.8,2)))
# 
# # SES
# ses=sapply(1:12,FUN=function(x) {
#       ses= (y[1,x] - mean(y[,x], na.rm=T)) /sd(y[,x], na.rm=T)
#       return(ses)   
#      } )  
# 
# plot(x, ses, type='h', lwd=10, lend=1, col=c('grey', 'black'), ann=F, axes=F, xlim=c(1,12), ylim=c(-4,13))
# abline(h=0)
#   
# axis(1, at=c(1.5,3.5,5.5,7.5,9.5,11.5), line=-0.5, labels = rep(c('-','+'),3), tcl=-0.3, lty=0, cex.axis=1.2)
# axis(1, at=c(2.5,6.5,10.5), line=0.5,lty=0,hadj=0.5,labels = c('All plots','Grasslands', 'Woodlands'),cex.axis=0.8)
# axis(1, at=c(0.5,4.5,8.5,12.5), lwd=0,lwd.ticks=1, labels = NA)
# axis(1, at=c(2.5,6.5,10.5), lwd=0,lwd.ticks=1, labels = NA, tcl=-0.2)
# axis(2, las=1, cex.axis=0.8, mgp=c(3,0.7,0), tcl=-0.4 )
# # title(ylab='Standardized frequency', line=2.5, cex.lab=0.9)
# abline(v=4.5, lty='dotted')                                                    
# abline(v=8.5, lty='dotted')
# box(bty='o')  
# sig=apply(cbind(sig.inf, sig.sup),1, min, na.rm=T)
# sig=sapply(sig, p2star)
# sig[sig=='ns']=''
# #mtext(side=3,line=-1, at=x, text=sig, cex.lab=0.7)    
# text(y=ses[ses>0]+0.5, x=x[ses>0], label=sig[ses>0], cex.lab=0.7)
# text(y=ses[ses<0]-0.5, x=x[ses<0], label=sig[ses<0], cex.lab=0.7)
# mtext(side=1,line=-1, at=1:12, text=y[1,], cex=0.6, col='grey40', font=3)
# mtext(side=3,line=0,adj=0.5,  text="b) Native richness", cex=0.9)
# 
# # ALIEN richness
# names=c("Natives.negSRali", "Aliens.negSRali","Natives.posSRali", "Aliens.posSRali")    
# y=data.frame(SR=shuff.all[1:1000,names], grass= shuff.grass[1:1000,names],wood=shuff.wood[1:1000,names])
#     # significant differences ?
#    sig.inf=sapply(1:12,FUN=function(x) {
#       p=(sum((y[,x]<y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
#       return(p)   
#      } )      
#                                          
#   sig.sup=sapply(1:12,FUN=function(x) {
#        p=(sum((y[,x]>y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
#       return( p )   
#      } )      
# 
# # SES
# ses=sapply(1:12,FUN=function(x) {
#       ses= (y[1,x] - mean(y[,x], na.rm=T)) /sd(y[,x], na.rm=T)
#       return(ses)   
#      } )  
# 
# plot(x, ses, type='h', lwd=10, lend=1, col=c('grey', 'black'), ann=F, axes=F, xlim=c(1,12), ylim=c(-4,13))
# abline(h=0)
#   
# axis(1, at=c(1.5,3.5,5.5,7.5,9.5,11.5), line=-0.5, labels = rep(c('-','+'),3), tcl=-0.3, lty=0, cex.axis=1.2)
# axis(1, at=c(2.5,6.5,10.5), line=0.5,lty=0,hadj=0.5,labels = c('All plots','Grasslands', 'Woodlands'),cex.axis=0.8)
# axis(1, at=c(0.5,4.5,8.5,12.5), lwd=0,lwd.ticks=1, labels = NA)
# axis(1, at=c(2.5,6.5,10.5), lwd=0,lwd.ticks=1, labels = NA, tcl=-0.2)
# axis(2, las=1, cex.axis=0.8, mgp=c(3,0.7,0), tcl=-0.4 )
# # title(ylab='Standardized frequency', line=2.5, cex.lab=0.9)
# abline(v=4.5, lty='dotted')                                                    
# abline(v=8.5, lty='dotted')
# box(bty='o')  
# sig=apply(cbind(sig.inf, sig.sup),1, min, na.rm=T)
# sig=sapply(sig, p2star)
# sig[sig=='ns']=''
# #mtext(side=3,line=-1, at=x, text=sig, cex.lab=0.7)    
# text(y=ses[ses>0]+0.5, x=x[ses>0], label=sig[ses>0], cex.lab=0.7)
# text(y=ses[ses<0]-0.5, x=x[ses<0], label=sig[ses<0], cex.lab=0.7)
# mtext(side=1,line=-1, at=1:12, text=y[1,], cex=0.7, col='grey40', font=3)
# mtext(side=3,line=0,adj=0.5,  text="c) Alien richness", cex=0.9)
# 
# 
# mtext(2, text='Standardized frequency', outer=T, line=1, cex=0.8, las=0)
# 
# 
# 
# # autre orga

 # By subcomponent
plot.sim2=function(y, space=0.9){
    # significant differences ?
       sig.inf=sapply(1:4,FUN=function(x) {
          p=(sum((y[,x]<y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
          return(p)   
         } )      
                                             
      sig.sup=sapply(1:4,FUN=function(x) {
           p=(sum((y[,x]>y[1,x]), na.rm=T) + sum(y[,x]==y[1,x], na.rm=T)/2)/dim(y)[1]
          return( p )   
         } )      
    # SES
    ses=sapply(1:4,FUN=function(x) {
          ses= (y[1,x] - mean(y[,x], na.rm=T)) /sd(y[,x], na.rm=T)
          return(ses)   
         } )  
    
    plot(x, ses, type='h', lwd=18, lend=1, col=c('grey', 'black'), ann=F, axes=F, xlim=c(1,4),
         ylim=ylim)
    abline(h=0)
    axis(1, at=c(1.75,3.25), line=-1.2,adj=0.5, labels = c('negative','positive'), tcl=-0.3, lty=0, cex.axis=0.7)
#   axis(1, at=c(2.5), lwd=0,lwd.ticks=1, labels = NA, tcl=0.2)
    abline(v=2.5, lty='dotted', col="grey")                                                    
#   box(bty='l')
    
    sig=apply(cbind(sig.inf, sig.sup),1, min, na.rm=T)
    sig=sapply(sig, p2star)
    sig[sig=='ns']=''
    text(y=ses[ses>0]+space, x=x[ses>0], label=sig[ses>0], cex.lab=0.7)
    text(y=ses[ses<0]-space, x=x[ses<0], label=sig[ses<0], cex.lab=0.7)
    mtext(side=1,line=-1, at=x, text=y[1,], cex=0.7, col=c('grey40'), font=0)
    axis(2, las=1, cex.axis=0.7, mgp=c(3,0.2,0), tcl=0.1 )
}

# graph
par(mfrow=c(3,3), oma=c(1,6,2,1), mar=c(2,1,0,0), cex=0.9, las=1)
x=c(1.5,2,3,3.5)
ylim=c(-5,6)

# total richness
names=c("Natives.negSR","Aliens.negSR","Natives.posSR","Aliens.posSR")
plot.sim2(y=shuff.all[1:1000,names])
# axis(2, las=1, cex.axis=0.65, mgp=c(3,0.5,0), tcl=- 0.3)
plot.sim2(y=shuff.grass[1:1000,names])
plot.sim2(y=shuff.wood[1:1000,names])
legend(x=2.6, y=5.5, legend=c('Natives','Aliens'), fill=c('grey', 'black'),border=NA, cex=0.8, bty='n')

# Native richness
names=c("Natives.negSRnat", "Aliens.negSRnat", "Natives.posSRnat" , "Aliens.posSRnat")
plot.sim2(y=shuff.all[1:1000,names])
# axis(2, las=1, cex.axis=0.65, mgp=c(3,0.5,0), tcl=- 0.3)
plot.sim2(y=shuff.grass[1:1000,names])
ylim=c(-5,6.5)
plot.sim2(y=shuff.wood[1:1000,names])

       # ALIEN richness
names=c("Natives.negSRali", "Aliens.negSRali","Natives.posSRali", "Aliens.posSRali") 
ylim=c(-5,6)
plot.sim2(y=shuff.all[1:1000,names])
# axis(2, las=1, cex.axis=0.65, mgp=c(3,0.5,0), tcl=- 0.3)
plot.sim2(y=shuff.grass[1:1000,names])
ylim=c(-5.8,15)
plot.sim2(y=shuff.wood[1:1000,names], space=1.4)

mtext(2, at=0.54,text='Standardized frequency', outer=T, line=0.2, cex=0.8, las=0)
mtext(side=3, at=c(0.18,0.52,0.85), text=c('All plots', 'Grasslands','Woodlands'),
      outer=T, las=0, line=0.5, cex=0.9)    
mtext(side=2, at=c(0.21,0.54,0.88), text=c('Alien\nrichness', 'Native\nrichness', 'Total\nrichness'),
      outer=T, las=1, line=4, adj=0.5, font=0, cex=0.9)
mtext(side=2, at=c(0.15,0.48,0.82), text=c('(c)', '(b)', '(a)'),
      outer=T, las=1, line=4, adj=0.54, font=0, cex=0.8)
 