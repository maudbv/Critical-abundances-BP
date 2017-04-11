# Phylogeny of BP species
library(ape, picante)
 load("saved Rdata/lmaphy.Rdata")
# load("saved Rdata/Banks peninsula dataset.Rdata")
# load("saved Rdata/lmaphy.Rdata")
# load("saved Rdata/article 1/effect summary.Rdata")



##  adding the missing tips

misstip=unique(species[!species$tip%in%lmaphy$tip.label ,c("Sp.code","Genus","Species","tip")]) # 456 missing species

length(which(misstip$Genus %in% lmaphy$node.label)) #456 unplaced species

modif.lmaphy=lmaphy
# zoom(modif.lmaphy, c(40:60, 5471), show.node.label=T)
# l=subtrees(modif.lmaphy, wait=T)
# which(modif.lmaphy$node=="Eupolypods")
# plot(l[[26]])

for (i in 1:length(misstip$tip)) {
if  (misstip$Genus[i]%in%modif.lmaphy$node.label) {
   l=mean(unique(modif.lmaphy$edge.length[which(modif.lmaphy$edge[,1]== 
     (which(modif.lmaphy$node==misstip$Genus[i])+length(modif.lmaphy$tip.label)))]))
   t=paste("(",misstip$tip[i],":",l,");",sep="")
   t=read.tree(text=t)
   t$root.edge=1
   modif.lmaphy=bind.tree(modif.lmaphy,t,where=which(modif.lmaphy$node==misstip$Genus[i])+length(modif.lmaphy$tip.label),position=0)
   }
 } 



misstip=species[!species$tip%in%modif.lmaphy$tip.label ,] # now only 197 species unplaced
missfam=unique(as.character(misstip[which(!misstip$Family%in%modif.lmaphy$node),]$Family)) # 14 unplaced native families only
unique(misstip[which(!misstip$Family%in%modif.lmaphy$node),c("Genus","Family")]) # 17 missing genera

# missing 6 of the "important" species with significant effects = 5 natives and 1 alien
# but only 1 native with negative effect on total richness
# And only 1 native with an unplaced family...
missimp=na.omit(effect.summary[effect.summary$SpeciesCode %in% misstip$Sp.code & effect.summary$total.effect>0,])

## best tree for BP :
tmp=modif.lmaphy$tip.label[!modif.lmaphy$tip.label %in%species$tip]
BPtree=drop.tip(modif.lmaphy,tmp,trim.internal=T)

# add subspecies info
tmp=which(duplicated(species$tip))

for (i in tmp) {
  if  (species$tip[i]%in%BPtree$tip.label) {
    l=mean(unique(BPtree$edge.length[which(BPtree$edge[,1]== 
                                                   (which(BPtree$node==species$Genus[i])+length(BPtree$tip.label)))]))
    t=paste("(",species$tipcomplete[i],":",l,");",sep="")
    t=read.tree(text=t)
    t$root.edge=1
    BPtree=bind.tree(BPtree,t,where=which(BPtree$node==species$Genus[i])+length(BPtree$tip.label),position=0)
  }
} 






BPtree=reorder(BPtree, order = "cladewise")

unique(species[!species$tip%in%modif.lmaphy$tip.label ,c("Sp.code","Genus","Species","tip")]) # 352 missing species




#save(lmaphy, BPtree, file="saved Rdata/phylogenies.Rdata")

# graphic
alientips=BPtree$tip.label[BPtree$tip.label%in% species$tip[species$ALIEN==1]]
x11()
plot.phylo(BPtree, type="fan", show.node.label=T, tip.color=c("black","indianred3")[as.numeric(BPtree$tip.label %in% alientips) +1 ],
           cex=0.45, no.margin=T)

par(mar=c(1,1,1,1))
plot(BPtree)
zoom(BPtree, focus=1:4, show.node.label=T,cex=0.8) 

plot(BPtree)
zoom(BPtree, focus=80:90, show.node.label=T,cex=0.8) 
