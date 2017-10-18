library(ROCR)


## ----- misc -----------
setwd('~/b_affy')
# setwd('F:\\Project_group\\Git_Project\\FHIR\\b_affy\\bn_mod')

library(bnlearn)
library(Rgraphviz)

data("coronary")


# structure learning
struc_learn <- function(dataset, algorithm='mmhc'){
  if(algorithm == 'mmhc'){
    net_stru = mmhc(dataset)
  } else if(algorithm == 'rsmax2'){
    net_stru = rsmax2(dataset)
  } else if(algorithm == 'tabu'){
    net_stru = tabu(dataset)
  } else if(algorithm == 'hc'){
    net_stru = hc(dataset)
  } else{
    stop(("algorithm name error !!"))
  }
  return(net_stru)
}

net_info <- function(net_stru){
  print(net_stru)
  bnlearn::graphviz.plot(net_stru)
}

# parameter learning 
param_learn <- function(net_stru, dataset){
  cpt <- bn.fit(net_stru, dataset)
  return(cpt)
}


### ----- operations ---------

# network learning

net_stru <- struc_learn(coronary, algorithm='rsmax2')
cpts <- param_learn(net_stru, coronary)
net_info(net_stru)



bt_strength <- boot.strength(data = coronary, algorithm = 'rsmax2')
pred <- as.prediction(bt_strength, net_stru)
perf <- performance(pred, 'tpr', 'fpr')
plot(perf, colorize = T)

