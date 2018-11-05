subs <- comm[  ,names(comm) %in% natives]

tbl <- t(apply(subs[1:10,], 1, FUN = function(j) {
  return(table(factor(x= subs[j, ], levels = 1:7)))
}))
