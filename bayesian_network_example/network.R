# setwd("F:\\Project_group\\Git_Project\\FHIR\\b_affy\\r_analyze")

library("bnlearn")
# library("BiocGenerics")

# structure learning
dataset = "test_data.txt"
datas = read.table(dataset, header = TRUE, colClasses = c("numeric"))
names(datas)<-c("L", "1", "2", "3", "4", "5", "6", "7", "8")
mmhc_net = mmhc(datas)
print(mmhc_net)
graphviz.plot(mmhc_net)

# mmhc_net = set.arc(mmhc_net, from = "B", to = "A")

# parameter learning 
params = bn.fit(mmhc_net, datas)