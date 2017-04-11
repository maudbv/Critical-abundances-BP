p2star <- function(x, marginal=F, quantiles=NULL) {
  # transforming p values into stars
 # x : a p value 
  s<-NA
  x <- as.numeric(x)
  if(!is.na(x)){
  if (!is.null(quantiles)) q<-quantiles
  else q<-c(0.1,0.05,0.01,0.001)
  if (x>q[2]) s<-'ns'
  if(marginal) if (x<q[1]) s<-'a'
  if (x<=q[2]) s<-'*'
  if (x<=q[3]) s<-'**'
  if (x<=q[4]) s<-'***'
}
  return(s)
} 
