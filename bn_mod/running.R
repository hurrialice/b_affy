setwd("F:\\Project_group\\Git_Project\\FHIR\\b_affy\\test_branch")

# data pre process, including feature select and discrezation
source("pre_process.R")


dataset <- load_data(data_file = "feature_table.rds")
subset <- feature_filter(dataset, correlation=0.8)
fealist <- feature_list(subset)
disset <- discrezation(subset, psigma=0.4)

# structure learning
source("network.R")
net_stru <- struc_learn(disset, algorithm='mmhc')
net_info(net_stru)
