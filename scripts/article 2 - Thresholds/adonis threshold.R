#### PERMANOVA  and CAPSCALE : 
# Beta-dissimilarity of invaded communities below and above threshold of target species


#### Species Matching above/below threshold for ALL SPECIES   ###########


dbrda.SM.th<- list()
adonis.SM.th <- list()
betadisper.SM.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,colSums(community)>0]
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  D <- 1-SMsim(community)
  
  dbrda.SM.th[[i]]  <- capscale(D ~ var)  
  adonis.SM.th[[i]] <- adonis(D ~ var)
  betadisper.SM.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.SM.th,adonis.SM.th,betadisper.SM.th, file = "saved Rdata/article 2 - threshold/dbRDA_SM.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_SM.Rdata")

x11()
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.SM.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  ymax= max(abs(scores(dbrda.results, display = "sites")[,2]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-ymax,ymax), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.SM.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.SM.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 

#### SMabove/below threshold for NATIVE SPECIES   ###########

dbrda.nat.SM.th<- list()
adonis.nat.SM.th <- list()
betadisper.nat.SM.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  D <- 1-SMsim(community)
  
  
  dbrda.nat.SM.th[[i]]  <- capscale(D ~ var)  
  adonis.nat.SM.th[[i]] <- adonis(D ~ var)
  betadisper.nat.SM.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.nat.SM.th,adonis.nat.SM.th,betadisper.nat.SM.th, file = "saved Rdata/article 2 - threshold/dbRDA_SM.natives.Rdata" )

#ordination plot
load("saveb d Rdata/article 2 - threshold/dbRDA_adj.bray.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.nat.SM.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xrange= range(scores(dbrda.results, display = "sites")[,1])
  yrange= range(scores(dbrda.results, display = "sites")[,2])
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=yrange, xlim=xrange) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.nat.SM.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.nat.SM.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 

#### SM above/below threshold for ALIEN SPECIES   ###########
## = using adjusted bray curtic distances with dummy variable


dbrda.ali.SM.th<- list()
adonis.ali.SM.th <- list()
betadisper.ali.SM.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% aliens] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  D <- 1-SMsim(community)
  
  dbrda.ali.SM.th[[i]]  <- rda(D ~ var)  
  adonis.ali.SM.th[[i]] <- adonis(D ~ var, method = "euclid")
  betadisper.ali.SM.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.ali.SM.th,adonis.ali.SM.th,betadisper.ali.SM.th, file = "saved Rdata/article 2 - threshold/dbRDA_SM.aliens.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_SM.aliens.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.ali.SM.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% aliens] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  ymax= max(abs(scores(dbrda.results, display = "sites")[,2]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-ymax,ymax), xlim=c(-xmax,xmax)) 
  
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.ali.SM.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.ali.SM.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 




#### With null models  :  ###############
#### Species Matching above/below threshold for ALL SPECIES   ###########


dbrda.SM.th<- list()
adonis.SM.th <- list()
betadisper.SM.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  D <- 1-SMsim(community)
  
  dbrda.SM.th[[i]]  <- rda(D ~ var)  
  adonis.SM.th[[i]] <- adonis(D ~ var)
  betadisper.SM.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.SM.th,adonis.SM.th,betadisper.SM.th, file = "saved Rdata/article 2 - threshold/dbRDA_SM.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_SM.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.SM.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  ymax= max(abs(scores(dbrda.results, display = "sites")[,2]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-ymax,ymax), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.SM.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.SM.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 

#### SM above/below threshold for NATIVE SPECIES   ###########

load( file="saved Rdata/article 2 - threshold/dist.SM_Natives_499reps.Rdata")

dbrda.nat.SM.th<- list()
adonis.nat.SM.th <- list()
betadisper.nat.SM.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
 
  D <- dist.SM.SRnat[[i]]
  
  
  dbrda.nat.SM.th[[i]]  <- capscale(D ~ var)  
  adonis.nat.SM.th[[i]] <- adonis(D ~ var)
  betadisper.nat.SM.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.nat.SM.th,adonis.nat.SM.th,betadisper.nat.SM.th, file = "saved Rdata/article 2 - threshold/dbRDA_SM.natives.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_SM.natives.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.nat.SM.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  ymax= max(abs(scores(dbrda.results, display = "sites")[,2]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-ymax,ymax), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.nat.SM.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.nat.SM.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 

#### SM above/below threshold for ALIEN SPECIES   ###########
## = using adjusted bray curtic distances with dummy variable


dbrda.ali.SM.th<- list()
adonis.ali.SM.th <- list()
betadisper.ali.SM.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% aliens] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  D <- 1-SMsim(community)
  
  dbrda.ali.SM.th[[i]]  <- rda(D ~ var)  
  adonis.ali.SM.th[[i]] <- adonis(D ~ var, method = "euclid")
  betadisper.ali.SM.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.ali.SM.th,adonis.ali.SM.th,betadisper.ali.SM.th, file = "saved Rdata/article 2 - threshold/dbRDA_SM.aliens.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_SM.aliens.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.ali.SM.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% aliens] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  ymax= max(abs(scores(dbrda.results, display = "sites")[,2]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-ymax,ymax), xlim=c(-xmax,xmax)) 
  
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.ali.SM.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.ali.SM.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 



#### Euclidean presence absence above/below threshold for ALL SPECIES   ###########


dbrda.E.th<- list()
adonis.E.th <- list()
betadisper.E.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
     
     
  dbrda.E.th[[i]]  <- rda(community ~ var)  
  adonis.E.th[[i]] <- adonis(community ~ var, method = "euclid")
  betadisper.E.th [[i]]<- betadisper(vegdist(community , method = "euclid"), group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.E.th,adonis.E.th,betadisper.E.th, file = "saved Rdata/article 2 - threshold/dbRDA_euclid.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_adj.bray.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.E.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-2,2), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.E.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.E.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 

#### Euclidean presence absence above/below threshold for NATIVE SPECIES   ###########


dbrda.nat.E.th<- list()
adonis.nat.E.th <- list()
betadisper.nat.E.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  #   
  #   D <- vegdist(community>0, method = "euclid")
  
  dbrda.nat.E.th[[i]]  <- rda(community ~ var)  
  adonis.nat.E.th[[i]] <- adonis(community ~ var, method = "euclid")
  betadisper.nat.E.th [[i]]<- betadisper(vegdist(community , method = "euclid"), group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.nat.E.th,adonis.nat.E.th,betadisper.nat.E.th, file = "saved Rdata/article 2 - threshold/dbRDA_euclid.natives.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_adj.bray.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.nat.E.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-2,2), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.nat.E.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.nat.E.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 

#### Euclidean presence absence above/below threshold for ALIEN SPECIES   ###########
## = using adjusted bray curtic distances with dummy variable


dbrda.ali.E.th<- list()
adonis.ali.E.th <- list()
betadisper.ali.E.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% aliens] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  #   
  #   D <- vegdist(community>0, method = "euclid")
  
  dbrda.ali.E.th[[i]]  <- rda(community ~ var)  
  adonis.ali.E.th[[i]] <- adonis(community ~ var, method = "euclid")
  betadisper.ali.E.th [[i]]<- betadisper(vegdist(community , method = "euclid"), group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.ali.E.th,adonis.ali.E.th,betadisper.ali.E.th, file = "saved Rdata/article 2 - threshold/dbRDA_euclid.aliens.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_adj.bray.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.ali.E.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% aliens] 
  community <- ceiling(community>0)
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-2,2), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.ali.E.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.ali.E.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 



### Raup-crick distances on presence - absence of NATIVES ##########
save(impsp, comm, natives, realgrasslands, myraupcrick, file ="data for raupcrick.Rdata")


setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
library(parallel)
library(doParallel)
cl <- makeCluster(6)
system.time(
  dist.RC.SRnat <- parLapplyLB( cl = cl, X = 1:11, fun = function(i) {
    load(file ="C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/data for raupcrick.Rdata")
    library(vegan)
    memory.limit(4095)
    sp=impsp[i]
    # select plots of grassland where the target is present
    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
    community <- community[,names(community) %in% natives]  # keep only natives
    community <- community[,which(colSums(community) >0) ]  # remove natives which never cooccur with sp
    
    dist <- myraupcrick(community, nreps=499)
    # dist <- raupcrick(community, nsimul=999) #it is a bit faster than mine but has pb with allocating memory
    gc()
    print(i)
    return(dist)
  })
)
stopCluster(cl)
save(dist.RC.SRnat, file="saved Rdata/article 2 - threshold/dist.RC_Natives_499reps.Rdata")


## Ordination with Raup crick

load( file="saved Rdata/article 2 - threshold/dist.RC_Natives_499reps.Rdata")
load(file ="C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/data for raupcrick.Rdata")
load("saved Rdata/article 2 - threshold/article threshold 1.2.Rdata")

dbrda.rc.SRnat.th<- list()
adonis.rc.SRnat.th <- list()
betadisper.rc.SRnat.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives]  # keep only natives
  community <- community[,which(colSums(community) >0) ]  # remove natives which never cooccur with sp
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  # D <- dist.RC[[i]]
  
  bcom <- ceiling(community/10) 
  tri <- matrix(FALSE, dim(bcom)[1], dim(bcom)[1])
  tri <- row(tri) > col(tri) # select the lower triangle of the matrix, excluding diagonal values
  D <- matrix(NA,  dim(bcom)[1], dim(bcom)[1])
  D[tri] <-dist.RC.SRnat[[i]]
  D <- as.dist(D)
  
  dbrda.rc.SRnat.th[[i]]  <- capscale(D ~ var)  
  adonis.rc.SRnat.th[[i]] <- adonis(D ~ var)
  betadisper.rc.SRnat.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}
save(dbrda.rc.SRnat.th,adonis.rc.SRnat.th,betadisper.rc.SRnat.th, dist.RC, file = "saved Rdata/article 2 - threshold/dbRDA_RC.SRnat.Rdata" )

#ordination plot
load(file = "saved Rdata/article 2 - threshold/dbRDA_RC.SRnat.Rdata" )


par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.rc.SRnat.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives]  # keep only natives
  community <- community[,which(colSums(community) >0) ]  # remove natives which never cooccur with sp
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  ymax= max(abs(scores(dbrda.results, display = "sites")[,2]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-ymax, ymax), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  #   
  #   ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
  #            groups=treat[ treat == "below"],draw="polygon",
  #            col="white",label=F)
  #   
  #   ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
  #            groups=treat[ treat == "above"],draw="polygon",
  #            col="grey30",label=F)
  
  ordiellipse(scores(dbrda.results, display = "sites")[grep("below",treat),],
              groups=treat[ treat == "below"],draw="polygon", kind="sd",conf = 0.95,
              col="white",label=F)
  
  ordiellipse(scores(dbrda.results, display = "sites")[grep("above",treat),],
              groups=treat[ treat == "above"],draw="polygon",kind="sd",conf = 0.95,
              col="grey30",label=F)
  
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.rc.SRnat.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.rc.SRnat.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 

### Raup-crick distances on presence - absence of ALIENS ##########
save(impsp, comm,aliens, natives, realgrasslands, myraupcrick, file ="data for raupcrick.Rdata")


setwd("C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/")
library(parallel)
library(doParallel)
cl <- makeCluster(6)
system.time(
  dist.RC.SRali <- parLapplyLB( cl = cl, X = 1:11, fun = function(i) {
    load(file ="C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/data for raupcrick.Rdata")
    library(vegan)
    memory.limit(4095)
    sp=impsp[i]
    community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
    community <- community[,names(community) %in% aliens] 
    
    dist <- myraupcrick(community, nreps=499)
    # dist <- raupcrick(community, nsimul=999) #it is a bit faster than mine but has pb with allocating memory
    gc()
    print(i)
    return(dist)
  })
)
stopCluster(cl)
save(dist.RC.SRali, file="saved Rdata/article 2 - threshold/dist.RC_Aliens_499reps.Rdata")


## Ordination with Raup crick

load( file="saved Rdata/article 2 - threshold/dist.RC_Aliens_499reps.Rdata")
load(file ="C:/Users/bernarm2/Dropbox/Work/doc boulot/post doc Lincoln/R/data for raupcrick.Rdata")
load("saved Rdata/article 2 - threshold/article threshold 1.2.Rdata")

dbrda.rc.SRali.th<- list()
adonis.rc.SRali.th <- list()
betadisper.rc.SRali.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% aliens]  # keep only aliens
  community <- community[,which(colSums(community) >0) ]  # remove natives which never cooccur with sp
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  # D <- dist.RC[[i]]
  
  bcom <- ceiling(community/10) 
  tri <- matrix(FALSE, dim(bcom)[1], dim(bcom)[1])
  tri <- row(tri) > col(tri) # select the lower triangle of the matrix, excluding diagonal values
  D <- matrix(NA,  dim(bcom)[1], dim(bcom)[1])
  D[tri] <-dist.RC.SRali[[i]]
  D <- as.dist(D)
  
  dbrda.rc.SRali.th[[i]]  <- capscale(D ~ var)  
  adonis.rc.SRali.th[[i]] <- adonis(D ~ var)
  betadisper.rc.SRali.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}
save(dbrda.rc.SRali.th,adonis.rc.SRali.th,betadisper.rc.SRali.th, dist.RC.SRali, file = "saved Rdata/article 2 - threshold/dbRDA_RC.SRali.Rdata" )

#ordination plot
load(file = "saved Rdata/article 2 - threshold/dbRDA_RC.SRali.Rdata" )


par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.rc.SRali.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,names(community) %in% aliens] 
  community <- community[,which(colSums(community) >0) ]  # remove natives which never cooccur with sp
  
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  ymax= max(abs(scores(dbrda.results, display = "sites")[,2]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-ymax, ymax), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.rc.SRali.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.rc.SRali.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 




###### dbRDA with Bray-Curtis above/below threshold = ALL SPECIES #########

dbrda.SR.th<- list()
adonis.SR.th <- list()
betadisper.SR.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which(rownames(comm) %in% realgrasslands),-grep(sp, names(comm))]
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  D <- vegdist(community, method = "bray")
  
  dbrda.SR.th[[i]]  <- capscale(D ~ var)  
  adonis.SR.th[[i]] <- adonis(D ~ var)
  betadisper.SR.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}
save(dbrda.SR.th,adonis.th,betadisper.th, file = "saved Rdata/article 2 - threshold/dbRDA_bray.Rdata" )

#ordination plot
load(file = "saved Rdata/article 2 - threshold/dbRDA_bray.Rdata" )

par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.SR.th[[i]]
  sp=impsp[i]
  community <- comm[which(rownames(comm) %in% realgrasslands),-grep(sp, names(comm))]
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  plot(dbrda.results, display = "sites", type="points", axes=F, ann=F) 
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  points(scores(dbrda.results, display = "sites")[grep("above",treat),], 
         pch =21, cex= 0.75, col ="grey30", bg ="grey30", ann=F) 
  
  #   hulls<- ordihull(scores(dbrda.results, display = "sites"),groups = treat,draw="polygon", label =F)
  #   
  #   summary(hulls)
  #   points(summary(hulls)[1:2,1:2], pch = 22,bg= c("black","white"))
  #   
  f<- adonis.SR.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.SR.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 

#### ADJUSTED BRAY CURTIS above/below threshold for NATIVE SPECIES   ###########
## = using adjusted bray curtic distances with dummy variable
dist.adjBC <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% natives] 
  
  # add dummy variable to comm
  community$dummy <- 1
  dist.adjBC[[i]] <- vegdist(community, method = "bray")
  print(i)
}


dbrda.th<- list()
adonis.th <- list()
betadisper.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,names(community) %in% natives] 
  community$dummy <- 1
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  D <- dist.adjBC[[i]]
  
  dbrda.th[[i]]  <- capscale(D ~ var)  
  adonis.th[[i]] <- adonis(D ~ var)
  betadisper.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}

save(dbrda.th,adonis.th,betadisper.th, file = "saved Rdata/article 2 - threshold/dbRDA_adj.bray.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_adj.bray.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,names(community) %in% natives] 
  community$dummy <- 1
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-2,2), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 

#### ADJUSTED BRAY CURTIS above/below threshold for ALIEN SPECIES   ###########
## = using adjusted bray curtic distances with dummy variable
dist.SRali.adjBC <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),-grep(sp, names(comm))]
  community <- community[,names(community) %in% aliens] 
  
  # add dummy variable to comm
  community$dummy <- 1
  dist.adjBC[[i]] <- vegdist(community, method = "bray")
  print(i)
}


dbrda.SRali.th<- list()
adonis.SRali.th <- list()
betadisper.SRali.th <- list()
for (i in 1:length(impsp)) {
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,names(community) %in% aliens] 
  community$dummy <- 1
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  D <- dist.adjBC[[i]]
  
  dbrda.SRali.th[[i]]  <- capscale(D ~ var)  
  adonis.SRali.th[[i]] <- adonis(D ~ var)
  betadisper.SRali.th [[i]]<- betadisper(D, group = var,type ="centroid")
  print(paste(i, ":", sp))
}
save(dbrda.SRali.th,adonis.SRali.th,betadisper.SRali.th, file = "saved Rdata/article 2 - threshold/dbRDA_adj.bray.SRali.Rdata" )

#ordination plot
load("saved Rdata/article 2 - threshold/dbRDA_adj.bray.SRali.Rdata")
par(mfrow=c(3,4), mar = c(1,1,1,1))
for (i in 1:length(impsp)) {
  dbrda.results <- dbrda.SRali.th[[i]]
  sp=impsp[i]
  community <- comm[which((rownames(comm) %in% realgrasslands) & (comm[,sp]>0) ),]
  community <- community[,names(community) %in% aliens] 
  community$dummy <- 1
  
  var <- comm[rownames(community),sp]
  var[var < glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 0
  var[var >= glmSRnat.overall$impact.spread [sp,"th.CI"]] <- 1
  
  treat <- c("below", "above") [var+1]
  
  #   plot(dbrda.results, display = "sites", type="points", axes=T, ann=F) 
  
  xmax= max(abs(scores(dbrda.results, display = "sites")[,1]))
  plot(scores(dbrda.results, display = "sites"), 
       type="p", axes=F, ann=F, ylim=c(-2,2), xlim=c(-xmax,xmax)) 
  abline(v=0,h=0, lty="dotted")
  
  
  ordihull(scores(dbrda.results, display = "sites")[grep("below",treat),],
           groups=treat[ treat == "below"],draw="polygon",
           col="white",label=F)
  
  ordihull(scores(dbrda.results, display = "sites")[grep("above",treat),],
           groups=treat[ treat == "above"],draw="polygon",
           col="grey30",label=F)
  
  #   h2 <- ordispider(scores(dbrda.results, display = "sites"),
  #                  groups= var, spiders = c("centroid"),label=F)
  #   
  #points(x = unique(h2)[1] ,y = unique(h2)[2], pch = 22 )
  
  f<- adonis.SRali.th [[i]]
  mtext(3, text =paste("r2 = ", round(f$aov.tab$R2[1], 2),p2star(f$aov.tab$Pr[1] ), sep=""),
        line = -1, adj = 1, cex = 0.7)
  
  disp <- betadisper.SRali.th [[i]]
  mtext(3, text =paste("Disp",p2star(anova(disp)$Pr[1]), sep=" "),
        line = -2, adj = 1, cex = 0.7)
  
  mtext(3, text =sub("_", " ",species[sp, "tip"]), line= 0, font=3, cex = 0.8)
}

plot.new()
legend("center",legend = c(">critical abundance", "<critical abundance"), fill= c("grey30","white"), bty="n") 


