# build a [classifier]
setwd('/home/qingzhang/b_affy2/bn_mod')


library(readr)
library(bnlearn)
library(ROCR)
ft <- read_rds('featrue_table.rds')
pathways <- read_csv('raw/reactome_pathway.csv')

## ----- assist functions -----------
# reduce features

feature_filter_by_pathways <- function(ft, pa = pathways, pval.cut = 0.05, fdr.cut = 0.5, ngenes = 400){
  pathways <- pa[pa$`Entities pValue` < pval.cut & pa$`Entities FDR` < fdr.cut,]
  genes_in_pathway <- unname(sapply(pathways$`Submitted entities found`, function(vec) strsplit(vec, split = ';' )))

  genes_sorted <- names(sort(table(unlist(genes_in_pathway))))
  
  reduced_genes <- genes_sorted[1:min(ngenes, length(genes_sorted))]
  ft <- ft[, c(reduced_genes, 'label') ]
}



## make discrete by zscore, apply to single featrue table
discrezation <- function(dataset, psigma){

    label <- dataset[, which(colnames(dataset) == 'label')]
    dataset <- dataset[, -which(colnames(dataset) == 'label')] # remove 'label'
    
    out <- apply(dataset, 2, function(fvec){
      sigma <- sd(fvec)
      mean <- mean(fvec)
      fvec_scaled <- (fvec - mean)/sigma
      # assign -1, 0, 1
      fout_1 <- ifelse(abs(fvec_scaled) <= psigma, 0, 1)
      fout_2 <- ifelse(fvec_scaled < -psigma, -1, 1)
      apply(rbind(fout_1, fout_2), 2, min)
    })

    out <- cbind(out, label)
    out <- apply(out, 2, as.factor)
    as.data.frame(out)
  }



# define train and test set for subsequent five-fold cross validation
diff_train_from_test <- function(ft, k){
  # make sure all test/train data is a mixture of AMI and control
  label <- ft$label
  AMI_ids <- which(label == 1)
  
  
  AMI_ft <- ft[ft$label == 1,]
  control_ft <- ft[ft$label == 0,]
  
  ncon2test <- floor(21/k)
  nami2test <- floor(52/k) - ncon2test
  
  folds_data <- list()
  for (i in 1:k){
    control_test <- control_ft[seq(from = ncon2test*(i-1)+1, to = ncon2test*i),]
    ami_test <- AMI_ft[seq(from = nami2test*(i-1)+1, to = nami2test*i),]
    
    test <- rbind(ami_test, control_test)
    train <- ft[!(rownames(ft) %in% rownames(test)), ]
    
    folds_data[[i]] <- list(test, train)
    names(folds_data[[i]]) <- c('test', 'train')
  }
  
  folds_data # a nested list
}



# perform k-fold validation
# ft is discrete dataframe!
k_fold_validation <- function(ft, k){
  
  data <- diff_train_from_test(ft, k )
  #browser()
  auc <- vector()
  
  
  for (i in seq_along(data)){
    # define test from train
    train <- data[[i]]$train
    test <- data[[i]]$test
    
    # build model
    tan_stru <- tree.bayes(x = train, training = "label", explanatory = head(colnames(train), -1))
    fitted = bn.fit(tan_stru, train, method = "bayes")
    
    #browser()
    # make prediction to train data
    pred = predict(fitted, test, prob = TRUE)
    results_prob = data.frame(t(attributes(pred)$prob))
    #browser()
    # plot roc
    compare_pred <- prediction(results_prob$X1, test$label)
    print(table(pred, test$label))
    
    roc <- performance(compare_pred, 'tpr', 'fpr')
    png(filename = sprintf("roc_plots_30/roc_folds%d.png", i))
    ROCR::plot(roc, colorize=TRUE)
    dev.off()
    #browser()
    auc_misc <- performance(compare_pred, 'auc')
    auc[i] <- auc_misc@y.values[[1]]
    # standard cross-validation
    
  }
  cat('auc =', mean(auc))
  
}


# --------- operations ---------------
# filter featrues by their frequency in enriched pathways
ft <- feature_filter_by_pathways(ft, pa = pathways, pval.cut = 0.05, fdr.cut = 0.5, ngenes = 30)
# make feature table discrete
disc_ft <- discrezation(ft, psigma = 0.6)
# gather rplots for roc
dir.create('roc_plots_30/')
model_performance <- k_fold_validation(ft = disc_ft, k = 5)
