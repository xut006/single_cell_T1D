#selectLabel <- function(ctClust, selectBy){
#  
#}

ctSelect = ctClust
selectBy = c("RAD", "TNET", "PD1+/ICOS+")

for (i in selectBy){
  ctSelect <- ctSelect[ctSelect$age == i]
  
}

