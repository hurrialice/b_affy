setwd("F:\\Project_group\\Git_Project\\FHIR\\b_affy\\test_branch")

library("bnlearn")

# data pre process, including feature select and discrezation
source("pre_process.R")
data_file = "feature_table.rds"
dataset <- load_data(data_file)
subset <- feature_filter(dataset)
disset <- discrezation(subset)

# structure learning
struc_learn <- function(dataset, algorithm='mmhc'){
    if(algorithm == 'mmhc'){
        net_stru = mmhc(dataset)
    }
    return(net_stru)
}

net_info <- function(net_stru){
    print(net_stru)
    graphviz.plot(net_stru)
}

# parameter learning 
param_learn <- function(net_stru, dataset){
    cpt <- bn.fit(net_stru, dataset)
    return(cpt)
}