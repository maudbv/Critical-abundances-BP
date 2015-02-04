## trying a random forest method

library(party)

mydata=data.frame(comm[envplot$PLOTID,], SR=envplot$SR)


data.controls <- cforest_unbiased(ntree=1000, mtry=3)
set.seed(47)
data.ctree <- ctree(as.formula(paste("SR ~ " ,paste(aliens[1:4], collapse="+"))), data = mydata,controls=data.controls)
plot(data.ctree)
data.cforest <- cforest(as.formula(paste("SR ~ " ,paste(aliens[1:50], collapse="+"))), data = mydata,controls=data.controls)
data.cforest.varimp <- varimp(data.cforest, conditional =TRUE)

mydata=envplot[, c("SR","PANN0080","AVTEMP0080", "MINTEM0080","MAXTEM0080", "GDD",
                   "SOLRADR","MOISTR","UR1991_DNS","BLDG_DIST","DIST_HWSLD","DIST_METAL","DIST_PRIV",
                   "DIST_SRIV", "DEM_10", "ASPECT","SLOPE","PH_MID")]
data.controls <- cforest_unbiased(ntree=1000, mtry=5)
set.seed(47)
data.ctree <- ctree(SR ~. , data = mydata,controls=data.controls)

data.cforest <- cforest(SR ~. , data = mydata,controls=data.controls)
data.cforest.varimp <- varimp(data.cforest, conditional =TRUE)
dotplot(sort(data.cforest.varimp))

plot(data.ctree)
