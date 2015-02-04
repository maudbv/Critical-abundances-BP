## Test class-SR relationships against null models
nreps= 999

######## NULL MODEL 1 :  shuffling abundances within communities

sim1=simul1.class(db=databp,var="SR",min=10, alpha= 0.01, nreps=nreps); save(sim1, file='saved Rdata/article 1/null model results/sim1.Rdata')
sim1.nat=simul1.class(db=databp,var="SRnat", min=10,alpha=0.01, nreps=nreps); save(sim1.nat, file='saved Rdata/article 1/null model results/sim1.nat.Rdata')
sim1.ali=simul1.class(db=databp,var="SRali",min=10, alpha= 0.01, nreps=nreps); save(sim1.ali, file='saved Rdata/article 1/null model results/sim1.ali.Rdata')

sim1.grass=simul1.class(db=databp[databp$PlotName%in% grasslands,],var="SR",min=10, alpha= 0.01, nreps=nreps); save(sim1.grass, file='saved Rdata/article 1/null model results/sim1.grass.Rdata')
sim1.wood=simul1.class(db=databp[databp$PlotName%in% woodlands,],var="SR",min=10, alpha= 0.01, nreps=nreps); save(sim1.wood, file='saved Rdata/article 1/null model results/sim1.wood.Rdata')
sim1.ali.grass=simul1.class(db=databp[databp$PlotName%in% grasslands,],var="SRali",min=10, alpha= 0.01, nreps=nreps); save(sim1.ali.grass, file='saved Rdata/article 1/null model results/sim1.ali.grass.Rdata')
sim1.ali.wood=simul1.class(db=databp[databp$PlotName%in% woodlands,],var="SRali",min=10, alpha= 0.01, nreps=nreps); save(sim1.ali.wood, file='saved Rdata/article 1/null model results/sim1.ali.wood.Rdata')
sim1.nat.grass=simul1.class(db=databp[databp$PlotName%in% grasslands,],var="SRnati",min=10, alpha= 0.01, nreps=nreps); save(sim1.nat.grass, file='saved Rdata/article 1/null model results/sim1.nat.grass.Rdata')
sim1.nat.wood=simul1.class(db=databp[databp$PlotName%in% woodlands,],var="SRnat",min=10, alpha= 0.01, nreps=nreps); save(sim1.nat.wood, file='saved Rdata/article 1/null model results/sim1.nat.wood.Rdata')

######## SIMULATION 1 : rho values only !!

sim1.rhos=simul1.class.rhos(db=databp,var="SR",min=10, alpha= 0.01, nreps=nreps); save(sim1.rhos, file='sim1.rhos.Rdata')
sim1.rhos.nat=simul1.class.rhos(db=databp,var="SRnat",min=10, alpha= 0.01, nreps=nreps); save(sim1.rhos.nat, file='sim1.rhos.nat.Rdata')
sim1.rhos.ali=simul1.class.rhos(db=databp,var="SRali",min=10, alpha= 0.01, nreps=nreps); save(sim1.rhos.ali, file='sim1.rhos.ali.Rdata')

sim1.rhos.grass=simul1.class.rhos(db=databp[databp$PlotName%in% grasslands,],var="SR",min=10, alpha= 0.01, nreps=nreps); save(sim1.rhos.grass, file='sim1.rhos.grass.Rdata')
sim1.rhos.ali.grass=simul1.class.rhos(db=databp[databp$PlotName%in% grasslands,],var="SRali",min=10, alpha= 0.01, nreps=nreps); save(sim1.rhos.ali.grass, file='sim1.rhos.ali.grass.Rdata')
sim1.rhos.nat.grass=simul1.class.rhos(db=databp[databp$PlotName%in% grasslands,],var="SRnat",min=10, alpha= 0.01, nreps=nreps); save(sim1.rhos.nat.grass, file='sim1.rhos.nat.grass.Rdata')

sim1.rhos.wood=simul1.class.rhos(db=databp[databp$PlotName%in% woodlands,],var="SR",min=10, alpha= 0.01, nreps=nreps); save(sim1.rhos.wood, file='sim1.rhos.wood.Rdata')
sim1.rhos.ali.wood=simul1.class.rhos(db=databp[databp$PlotName%in% woodlands,],var="SRali",min=10, alpha= 0.01, nreps=nreps); save(sim1.rhos.ali.wood, file='sim1.rhos.ali.wood.Rdata')
sim1.rhos.nat.wood=simul1.class.rhos(db=databp[databp$PlotName%in% woodlands,],var="SRnat",min=10, alpha= 0.01, nreps=nreps); save(sim1.rhos.nat.wood, file='sim1.rhos.nat.wood.Rdata')


####### NULL MODEL 2 : random distribtuion of alien/native status

shuff.all=shuff.classSR(db = databp,species = species,plot=F, min.occur = 10,env=envplot,zeros=F, min.class = 3, alpha = 0.01, nreps = nreps) ; save(shuff.all, file='saved Rdata/article 1/null model results/shuff.all.Rdata')
shuff.grass=shuff.classSR(db = databp[databp$PlotName%in% grasslands,],plot=F,species = species,env=envplot,zeros=F,min.occur = 10, min.class = 3, alpha = 0.01, nreps = nreps) ; save(shuff.grass, file='saved Rdata/article 1/null model results/shuff.grass.Rdata')
shuff.wood=shuff.classSR(db = databp[databp$PlotName%in% woodlands,],species = species,plot=F,env=envplot,zeros=F, min.occur = 10, min.class = 3, alpha = 0.01, nreps = nreps) ; save(shuff.wood, file='saved Rdata/article 1/null model results/shuff.wood.Rdata')
