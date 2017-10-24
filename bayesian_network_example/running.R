setwd("F:/Project_group/Git_Project/FHIR/b_affy/test_branch")
setwd('/home/qingzhang/b_affy2/bayesian_network_example')

source("assistants_functions.R")

#library("readtext")
library("ROCR")
library("bnlearn")
library("Rgraphviz")
library(readr)

dataset <- load_data(data_file = "feature_table.rds")

net_string <- read_lines("net_structure.txt")     # net structure txt
net_stru <- model2network(net_string)      # set net structure by hand
# img <- graphviz.plot(net_stru, highlight = list(nodes = 'label')) # visualized net structure
sub_set <- subset(dataset, select = nodes(net_stru))    # build sub_set from net_stru
dis_set <- discretize(sub_set, method = "interval", breaks = 2) # discrezation by bnlearn 
train_set <- dis_set
test_set <- dis_set[21:35,]
cpts <- bn.fit(net_stru, train_set)   # fit net parameters, default method includes bayes and mle(default) 
predicts <- cpt_predict(test_set, cpts)
labels <- as.numeric(test_set$label)
pred <- prediction(predicts, labels)
roc <- performance(pred, 'tpr', 'fpr')
auc <- performance(pred, 'auc')


# sub_set <- sub_set_by_hand(dataset, start=1, end=5)    # build sub_set by appointed features
# dis_set <- discretize(sub_set, method = "interval", breaks = 2) # discrezation by bnlearn 
# net_stru <- struc_learn(sub_set)
# img <- graphviz.plot(net_stru, highlight = list(nodes = 'label')) # visualized net structure
# cpts <- bn.fit(net_stru, dis_set)   # fit net parameters, default method includes bayes and mle(default) 
# cpquery(cpts, event = (label == "(0.5,1]"), 
#               evidence = ((RABGAP1L  == "[7.73,8.86]") &
#                             (RPS11 == "(10.7,12]")))


# data("ROCR.simple")
# str(ROCR.simple)
# predicts <- ROCR.simple$predictions
# labels <- ROCR.simple$labels
# pred <- prediction(predicts, labels)
# roc <- performance(pred, 'tpr', 'fpr')
# auc <- performance(pred, 'auc')
plot(roc, colorize = T)
