# setwd("F:\\Project_group\\Git_Project\\FHIR\\b_affy\\test_branch")

library("bnlearn")
library("Rgraphviz")

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
    img <- graphviz.plot(net_stru, highlight = list(nodes = 'label'))
}


# parameter learning 
param_learn <- function(net_stru, dataset){
    cpt <- bn.fit(net_stru, dataset)
    return(cpt)
}

artificial_net <- function(net_string){
    return(model2network(net_string))
}