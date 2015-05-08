## palette couleur
couleurs=matrix(1:length(colours()),73,9)

### Identifier des jolies couleurs avec une fonction :
colorset=function()
{par(mar=c(2,2,0,0), cex.axis=0.6)
plot(1,1,type="n", xlim=c(1,10),ylim=c(1,73), ann=F, axes=F)
couleurs=matrix(1:length(colours()),73,9)
coord=NULL
for (i in 1:9)
{
  for (j in 1:73)
    {
    c=colours()[couleurs[j,i]]
    points(i,j,pch=22,col=c,bg=c)
    text (i,j,labels=c,pos=4, offset=0.3,cex = 0.6)
    text (i,j,labels=couleurs[j,i],pos=2, offset=0.3,cex = 0.6)
    coord=c(coord,i,j)
    }
}
colorset=identify(matrix(coord,length(coord)/2,2,byrow=T))
colorset=colours()[colorset]
return(colorset)
}

# exporter la liste des couleurs en jpeg
jpeg("figs/palette couleur.jpeg", width=1300, height=1200, pointsize=17 )
{
par(mar=c(2,2,0,0), cex.axis=0.6)
plot(1,1,type="n", xlim=c(1,10),ylim=c(1,73), ann=F, axes=F)
for (i in 1:9)
{
  for (j in 1:73)
    {
    c=colours()[couleurs[j,i]]
    points(i,j,pch=22,col=c,bg=c, cex=1.5)
    text (i,j,labels=c,pos=4, offset=0.3,cex = 0.6)
    text (i,j,labels=couleurs[j,i],pos=2, offset=0.3,cex = 0.6)
    }
}
}
dev.off()
