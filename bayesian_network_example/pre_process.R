# setwd("F:\\Project_group\\Git_Project\\FHIR\\b_affy\\test_branch")

# load data
load_data <- function(data_file="feature_table.rds", label="label", postive_label="AMI"){
    dataset <- readRDS(data_file)   # class data.frame
    # transform dataset["label"] from string label to num label 0 1
    for(i in 1:length(dataset[,label])){
        if(dataset[i,label] == postive_label){
            dataset[i,label] <- 1
        } else {
            dataset[i,label] <- 0
        }
    }
    dataset[,label] <- as.numeric(dataset[,label])
    # write.csv(dataset, "dataset.csv")
    return(dataset)
}

# feature selection
feature_filter <- function(dataset, label="label", fold_limit=1, correlation=0.8){
    # filter by fold change value
    # calculate correlation between sub_datasetfeature and outcome
    outcome <- dataset[,label]

    features <- list()
    idx <- 0
    for(i in 1:(length(dataset)-1)){
        f <- dataset[,i]
        cos <-(f %*% outcome)/(sqrt(f %*% f)*sqrt(outcome %*% outcome))
        if(abs(cos) > correlation){
            # cat("f: ", names(dataset)[i], "\ncos: ", cos)
            idx <- idx+1
            features[idx] <- names(dataset)[i]
        }
    }
    features[idx+1] <- label
    sub_set <- subset(dataset, select = unlist(features))
    # write.csv(subset, "subset.csv")
    return(sub_set)
}

sub_features <- function(dataset, start, end){
    if(start > end){
        tmp = end
        end = start
        start = tmp
    }
    outcome <- dataset$label
    subset <- dataset[,start:end]
    subset$label <- dataset$label
    return(subset) 
}

extract_features <- function(dataset, net_stru){
    return(subset(dataset, select = nodes(net_stru)))
}

feature_list <- function(dataset){
    fealist <- names(dataset)
    return(fealist)
}

# discrezation with thress threshold
discrezation <- function(dataset, label="label", psigma=0.3){

    # if(psigma > 0.5){
    #     psigma <- 1-psigma
    # }
    
    for(i in 1:length(dataset)){
        if(names(dataset)[i] == label){
            next
        }
        dataset[,i] <- dataset[,i] - mean(dataset[,i])
        for(j in 1:length(dataset[,i])){
            sigma <- var(dataset[,i])
            if(dataset[j,i] < (psigma*sigma)){
                dataset[j,i] <- -1
            } else {
                if(dataset[j,i] > (1-psigma)*sigma){
                    dataset[j,i] <- 1
                } else {
                    dataset[j,i] <- 0
                }
            }
        }
    }
    # write.csv(dataset, "dis_dataset.csv")
    return(dataset)
}