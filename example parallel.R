
#Parallel
#For On Local Computer
library(doParallel)
cl<-makeCluster(6)
registerDoParallel(cl)


output <-  foreach(i=1:1000, .combine='c') %dopar% {
  rnorm(1)
}

stopCluster(cl)

