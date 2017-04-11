### null model for threhold effects

library(vegan)

## shiflle abundance per site among observed species (SR, evenness, sp pool constant, only abundance of sp changes)
# a constrained shuffle of rows, with column sums/distribution not kept. 


fnm3=function(x){
xn= apply(x, 1, FUN=function(z){
   z[z>0]=sample(z[z>0])
   return(z)
 })
  return(as.matrix(t(xn))       
} ## continues to give an unexplained error message, but seems to work as it should...

fnm3=function(x=databp){
  xn=orderBy(~ PlotName + SpeciesCode, data=x)
 y=sapply(unique(x$PlotName), FUN=function(i) {
    z=xn[xn$PlotName==i,"abun"]
    xn[xn$PlotName==i,"abun"]<- sample(z, length(z))
  })
  
xn$abun=unlist(y)    
  return(xn)
} ## continues to give an unexplained error message, but seems to work as it should...


nm3=lapply(1:9, FUN=function(i) fnm3(x))

nm3.effect=lapply(nm3, FUN=function(x) {
  y=pairedclass.test(db=x)
  return(y$threshold)
})

null.thr= lapply(nm3.effect, FUN=function(x) {
  table(x$threshold)})