### Nomenclature matching 

resolve.names <- function(df=species, colname="tip",db=c("tNRS"), proxy=T,username=NULL, password=NULL){
require(httr)
library(taxize)

x= as.vector(df[,colname])
x=sub("_"," ", x)

# correct encodings:

Encoding(x)

x=sub("×" , "", x) #not an x but a multiply character in latin-1

## checking the inputs
if(!is.character(x)) warning("the name vector is not of character class")

proxy=readline(prompt = "use proxy? (y/n) ")
if (proxy=="y") {
##### Setting the PROXY configuration to query the online databases
if (is.null(username)) username <- readline(prompt = "Enter username for proxy: ")
if (is.null(password)) password <- readline(prompt = "Enter password for proxy: ")

reset_config() # erase previous configurations
set_config(
  use_proxy('proxy1.lincoln.ac.nz', 8080, username=username,password=password)
)
}s

# GNR resolve : only 349 at a time!!
# iplant
out=data.frame()
times=(floor(length(x)/359))
for( i in 0:times) {
 indices= (i*349+1) : min((i*349+349), length(x))
 y <- resolve(x[indices], db="iplant")
 out <- rbind(out,y )
}

# TNRS matching : Match taxonomic names using the Taxonomic Name Resolution Service (TNRS) 
# Gbif
if ("tnrs" %in% db) {
  tnrs.out<-  tnrs(query = as.vector(x),getpost="GET", source = "GBIF")
  tnrs.ncbi$submittedname = sub("\r", "",tnrs.ncbi$submittedname) # NCBI mysteriously adds "\r" suffixes
}

## taxonstand package for TPL
# tplck

## NZOR


corrected <- data.frame(
  originalname=df[,colname],
  submittedname = x,
  ncbi=tnrs.ncbi[match(as.vector(df[,colname]), tnrs.ncbi$submittedname), c("acceptedname", "score")],
  iplant=tnrs.iPlant[match(x, tnrs.iPlant$submittedname), c("acceptedname", "score")]
)

corrected$preferredname  <-  corrected$ncbi.acceptedname
corrected$preferredname[which(as.numeric(corrected$iplant.score)>0.9  & is.na(corrected$ncbi.score))] <-
corrected$iplant.acceptedname[which(as.numeric(corrected$iplant.score)>0.9 & is.na(corrected$ncbi.score))]

# reset configuration to normal
# if (proxy=="y") reset_config()

return(corrected)
}
