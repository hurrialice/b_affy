setwd("F:/Project_group/Git_Project/FHIR/b_affy/test_branch")

# data pre process, including feature select and discrezation
source("pre_process.R")
# structure learning
# struc_learn algorithm include hc, tabu, mmhc, rsmax2
source("network.R")

ana_by_fea <- function(dataset, img_root = "result_img_by_features/"){
    limit <- length(dataset)
    start <- 1
    end <- 1

    dir.create(img_root)

    while(end < limit){

        end <- start + 14
        if(end > limit){
            end <- limit
        }

        subset <- sub_features(dataset, start, end)
        dir_name <- paste(img_root,end,sep="")
        dir.create(dir_name)
        for (i in 1:100) {
            disset <- discrezation(subset, psigma=i/100)
            net_stru <- struc_learn(disset, algorithm='rsmax2')
            img <- net_info(net_stru)
            savePlot(paste(dir_name,"/",i,sep=""),"jpg")
        }
        start <- end
    }
}

sb_while <- function(dataset, i, dir_name){
    
    limit <- length(dataset)
    start <- 1
    end <- 1

        while(end < limit){

            end <- start + 14
            if(end > limit){
                end <- limit
            }

            subset <- sub_features(dataset, start, end)
            disset <- discrezation(subset, psigma=i/100)
            net_stru <- struc_learn(disset, algorithm='rsmax2')
            img <- net_info(net_stru)
            savePlot(paste(dir_name,"/",end,sep=""),"jpg")
            
            start <- end
        }
}

ana_by_sigma <- function(dataset, img_root = "result_img_by_sigma/"){
    

    dir.create(img_root)

    for (i in 1:100) {

        dir_name <- paste(img_root,i,sep="")
        dir.create(dir_name)
        sb_while(dataset, i, dir_name)
    }
}

library(readtext)

net_string <- readtext("net_structure.txt")
net_stru <- artificial_net(net_string$text)
img <- net_info(net_stru)

dataset <- load_data(data_file = "feature_table.rds")

sub_set <- extract_features(dataset, net_stru)
# subset <- sub_features(dataset, 1, 15) 
# # subset <- feature_filter(dataset, correlation=0.80)
# fealist <- feature_list(subset)

disset <- discrezation(sub_set, psigma=0.43)
# net_stru <- struc_learn(disset, algorithm='rsmax2')
# img <- net_info(net_stru)

library("ROCR")

cpts <- param_learn(net_stru, disset)
net_info(net_stru)

bt_strength <- boot.strength(data = disset, algorithm = 'rsmax2')
pred <- as.prediction(bt_strength, net_stru)
perf <- performance(pred, 'tpr', 'fpr')
auc <- performance(pred, 'auc')
plot(perf, colorize = T)


# ana_by_sigma(dataset)
# ana_by_fea(dataset)


