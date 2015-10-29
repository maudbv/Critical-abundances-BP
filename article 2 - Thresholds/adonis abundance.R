#### adonis on target abundance




####### ordination on species abundance #############

# abundance class colors
coltreat <-  colorRampPalette(c("palegoldenrod", "firebrick"))(7)
coltreat.tr <- paste(coltreat, "50", sep="")

dbrda.species<- list() adonis.species <- list() betadisper.species <- list()

for (i in 1:length(impsp)) { sp=impsp[i] community <-
comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp,
names(comm))] community <- community[,-which(colSums(community) == 0)] var <-
comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),sp]

dbrda.species[[i]]  <- capscale(community ~ var,  distance = "bray")
adonis.species[[i]] <- adonis(community ~ var, method = "bray")
betadisper.species [[i]]<- betadisper(vegdist(community, method = "bray"),
group = var,type ="centroid") print(paste(i, ":", sp)) }
save(dbrda.species,adonis.species,betadisper.species, file = "saved
Rdata/article 2 - threshold/dbRDA_bray.Rdata" )

#ordination plot
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.species[[i]]
  sp = impsp[i]

  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  treat =  comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),sp]
  plot(dbrda.results, display = "sites", type="points", axes=F, ann=F)

  for (j in sort(unique(treat))) {
    if (length(grep(j,treat)) >3) {
      #       huls <- ordihull(scores(dbrda.results, display = "sites")[grep(j,treat),],
      #        groups=treat[ treat == j],draw="polygon",
      #                        col=coltreat.tr[j+1],label=F)
      #
      ordiellipse(scores(dbrda.results, display = "sites")[grep(j,treat),],
                  groups=treat[ treat == j],draw="polygon",display="sites",
                  kind = c("sd"), conf=0.95,
                  col=coltreat.tr[j+1], label=F, border="grey30")
      print(j)

    }
  }

  #   ordihull(dbrda.results,groups=treat,draw="polygon", label =T)
  f<- adonis.species [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)

  disp <- betadisper.species [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)

  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}


plot(-1:8, -1:8, type ="n", ann=F, bty="n", axes=F)
for (i in 1:6) {
  polygon(x = c(i,i,i+1,i+1),   y =c(4,5,5, 4), col =coltreat[i], border =F)
}
text(1.5:6.6,rep(3.5,6),
     label = 1:6, cex = 0.7)

text(1.5:6.6,rep(3,6), srt=45, adj=1,
     label = c("Occasional" ,"Common-Occasional",  "Common",
               "Abundant-Common", "Abundant","Dominant"), cex = 0.7)

