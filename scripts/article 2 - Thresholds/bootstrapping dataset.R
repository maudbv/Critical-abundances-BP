### overall bootrapped data : Elena's advice


###-----------Creating bootstrap samples of the dataset :
bootstrap.dataset<- function(db=databp[databp$PlotName %in% realgrasslands,],
                             min.occur =5,  min.class = 2, nreps = 999) {
  
  ######### selecting species verifying conditions :
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=min.class)
                    &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
  db.modif <- db[which(db$SpeciesCode %in% a),]
  
  # list of species =to be targeted in analysis :
  sp.names <- unique(db.modif$SpeciesCode)
  
  
  ######### Bootstrapping the dataset : resampling plots with replacement
  
  # Unique set of plots to be resampled
  d = as.character(unique( db.modif$PlotName))
  
  # resampling the list of plot names with replacement *nreps* times
  boots <- sapply(1:nreps,FUN = function (k) sample(d, replace=TRUE))
  
  return(list(boots =boots,db.size=dim(db)[1],min.occur =min.occur,  min.class = min.class, nreps = nreps))
}

extract.indices <- function(boot.output, db = db) {
  
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=boot.output$min.occur)>=boot.output$min.class)
                    &  table(db$SpeciesCode, db$abun)[,1]>=boot.output$min.occur))
  db.modif <- db[which(db$SpeciesCode %in% a),]
  
  indices.table=NULL
  
  for (r in 1:boot.output$nreps){
    ## extract new datasets from original dataset using the bootstrapped line numbers:
    ind <-unlist(lapply(1:length(boot.output$boots[,r]), function (k)  which(db.modif$PlotName == boot.output$boots[k,r])))
    nrep = rep(r, length(ind))
    indices.table =rbind(indices.table, cbind(nrep, ind))
    print(paste(r, ":","(",Sys.time(),")"))
  }
  
  return(list(index = as.data.frame(indices.table), nreps = boot.output$nreps, min.class = boot.output$min.class, min.occur = boot.output$min.occur, db.size = boot.output$db.size))
}
