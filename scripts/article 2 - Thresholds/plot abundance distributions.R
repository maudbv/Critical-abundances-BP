M = glmSRnat.overall
db=databp[databp$PlotName %in% unimprovedgrasslands,]
variable = "SRnat"
  
  a <- row.names(M$dif)
  db.modif <- db[which(db$SpeciesCode %in% a),]
  db.modif$abun <- as.factor(  db.modif$abun)
  # list of species =to be targeted in analysis :
  sp.names <- a
 par( mfrow= c(9,8), mar = c(1,1,2,1))
 
 sp.names <- impsp
 par( mfrow= c(2,4), mar = c(3,2,2,1))
 
  for (i in 1:length(sp.names)) {
    
    sp <- sp.names[i]  # select species name
    print(paste(i, ":", sp, "(",Sys.time(),")"))
    
    col.line ="black"
    if (sp %in% aliens) col.line = "firebrick"
    

    ### FIRST : calculate GLM for observed dataset :
    sp.dat <- db.modif[as.character(db.modif$SpeciesCode)==sp,c("abun",variable,'PlotName')] # select occurrences of the species
    names(sp.dat) <- c('abun','var','PlotName')[1:dim(sp.dat)[2]]
    
    # plot(density(sp.dat$abun), ann = F, axes = F, xlim = c(1,6), col = col.line)
    barplot(table(sp.dat$abun), ann = F, xlim = c(0,7), ylim = c(0,250), col = col.line)
    
     box  (bty = "l")
    mtext(3, text=substr(paste(species[sp, "Genus"]," ", species[sp, "Species"], sep=""),1, 25), font = 3, outer= F,adj=0.1, cex=0.6, line=0.2, las = 1)
    
  }
    