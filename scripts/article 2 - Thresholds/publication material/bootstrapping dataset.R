### Two sequential functions performing the bootstrapping of plots across the vegetation survey dataset:

## DESCRIPTION ________________________________________________________

# These function uses a systematic vegetation survey of species abundances (categrotical abundance classes), and perfoms a resampling with replacement of all the communities/sites for the subset of sites where at least one of the focal species is present. 
# the identity of focal species is not given as an input, but selected based on criteria of frequency of occurrence and local abundances.
# The first funciton 'bootstrap.dataset' performs the resampling with replacement of all sites 'nreps' times, and provides a list of vectors of plot names
# The second function takes the output from the first function, and creates a dataframe of two columns: column one the number of the bootstrap loop (from 1 to nreps); column two contains the ordered resampled indices of plots for each bootstrap run. This output "bootindices" will be used to calculate the glms for each bootstrapped data set. 

## FUNCTION 1: Bootstrap dataset  ________________________________________________________

## INPUTS:
# db :        dataframe with the abundance of each species present in each plot
# min.occur:  minimum number of occurences of a focal species in a given abundance class
# min.class:  minimum number of classes a focal species needs to appear in sufficient occurrences (min.occur)
# nreps:      number of repetitions of bootstrapping

## OUTPUT is a list of the following dataframes:
# boots :     list of *nreps* ordered vectors of plot names, resampled with replacement from the original dataset
# db.size     size of the final dataset when selecting only sites where at least one focal species is present
# nb.plots    number of unique plots where at least one focal species is present
# focal.sp    vector of selected focal species (including both natives and alien species)
# min.occur   repeats the input
# min.class   repeats the input
# nreps       repeats the input


# FUNCTION 2: extract.indices  ________________________________________________________

## INPUTS:
# boots :     output from the function "bootstrap.dataset"
# db :        dataframe with the abundance of each species present in each plot

## OUTPUT is a list of the following dataframes:
# index:      two column dataframe: column 1 is the numebr of the bootstrap sample; column 2 is the ordered plot names
# db.size     size of the final dataset when selecting only sites where at least one focal species is present
# min.occur   repeats the input
# min.class   repeats the input
# nreps       repeats the input


# This function was written by Maud Bernard-Verdier - updated March 2018
# ______________________________________________________________________________________________________


###-----------Creating bootstrap samples of the dataset :
bootstrap.dataset<- function(db= dataset,
                             min.occur =5,  min.class = 2, nreps = 999) {
  
  ######### selecting focal species verifying conditions :
  a <- names(which( (rowSums(table(db$SpeciesCode, db$abun)[,2:6]>=min.occur)>=min.class)
                    &  table(db$SpeciesCode, db$abun)[,1]>=min.occur))
  
  ## Selecting sites where at least one of the focal species is present:
  db.modif <- db[which(db$SpeciesCode %in% a),]
 
  ######### Bootstrapping the dataset : resampling plots with replacement
  
  # Unique set of plots to be resampled
  d = as.character(unique( db.modif$PlotName))
  
  # resampling the list of plot names with replacement *nreps* times
  boots <- sapply(1:nreps,FUN = function (k) sample(d, replace=TRUE))
  
  return(list(boots =boots, db.size=dim(db)[1], nb.plots= length(d) , nb.focal = a, min.occur =min.occur,  min.class = min.class, nreps = nreps))
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
