### add my trait data

mytraits <- read.csv(file="data/my traits summary.csv", stringsAsFactor=F)

traitdata$mySLA <- mytraits[match(traitdata$Sp.code, mytraits$Sp.code), "SLA"]
traitdata$myH <- mytraits[match(traitdata$Sp.code, mytraits$Sp.code), "H"]
traitdata$mySM <- mytraits[match(traitdata$Sp.code, mytraits$Sp.code), "SM"]


# summarize available data:
traitdata$SLA= apply(!is.na(cbind(traitdata$SLA_NG,
                                  traitdata$SLA_LEDA,
                                  traitdata$SLA_ordo,
                                  traitdata$mySLA)),1,
                     function(x) as.numeric(sum(x)>0))
traitdata$Hmax= apply(!is.na(cbind(traitdata$H_NG,
                                   traitdata$H_LEDA,
                                   traitdata$Hmax_ordo,
                                   traitdata$PlantHeight_ecotrait,
                                   traitdata$myH)),1,
                      function(x) as.numeric(sum(x)>0))
traitdata$SM= apply(!is.na(cbind( traitdata$SM_LEDA,
                                  traitdata$SM_ordo,
                                  traitdata$SM_ecotrait,
                                  traitdata$mySM)),1,
                    function(x) as.numeric(sum(x)>0))