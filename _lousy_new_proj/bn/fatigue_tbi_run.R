# load libraries, set seeds, define functions
source('fatigue_net_assist.R')


# enter manually
data_type <- "tbi"
# read in phenodata, data matrix
load(paste(data_type,"backup.RData", sep = "_"))

# rename variables
cvars <- c("data.matrix", "ph", "tab")
lapply(cvars, function(cvar){
        stored_var <- get(paste(cvar, data_type, sep = "_"))
        assign(cvar, stored_var, envir = .GlobalEnv)
})

# discretlize data matrix
ds_ft <- make_discre(data.matrix, 0.6)

# LOOCV
tids <- make_par()


## initiation of void lists
init_list(clist = c('tabugs', 'true_fts', 'train_labels','predicts','labels',
                    'test_labels', 'tab_results', 'tabugsf'))

# for each fold,
# calculate differential expression table from train set,
# and extract featrue table.
for (i in seq_along(tids)){
        
        # define train and test
        tid <- tids[[i]]
        train <- data.matrix[,tid]
        
        # differential expression
        train_label <- as.factor(ph$disease_state[tid])
        tab <- DE_pipe(train, train_label)
        
        # give labels to DE genes
        tab <- anno(t = tab)
        tab_results[[i]] <- tab
        
        # modify label
        train_labels[[i]] <- ifelse(train_label == ph$disease_state[1], 1, 0 )
        test_labels[[i]] <- ifelse(as.factor(ph$disease_state[-tid]) == ph$disease_state[1], 1, 0)
        
        # real feature table for tabu search (train + test)
        true_ft <- shrink_features(ori_dm = ds_ft, fea_names = tab$PROBEID)
        true_fts[[i]] <- true_ft
        
        # boot tabu
        train <- true_ft[tid,] # note here, 'train' are discrete values
        tabugsf[[i]] <- future({boot_tabu(ori_tab = train)})
        
        
        print(paste0(i, 'th done at ', Sys.time()))
}

# evluate futures, could be time consuming...
tabug_boots <- lapply(tabugsf, FUN = value)
save(tabug_boots, tids, true_fts, file = paste(data_type, 'store.RData', sep = "_"))
paste0('bagging finished at ', Sys.time())

##### lazy method ##### 
# if you do not wish to wait for bootstraping, 
# just load this (result for one fold)...

# note that these data is not the same as the generated one 
# since I used different seed for fold generalization

# save(tabug_boots, tids, true_fts, file = 'lazy_load_test_data.RData')

#######################



## ------- model assessment  -----------
for (m in seq_along(tabug_boots)){
        
        # re-init fold info
        tabu_boot <- tabug_boots[[m]]
        tabug <- averaged.network(tabu_boot, threshold = 0.7) # ad hoc value
        true_ft <- true_fts[[m]]
        tid <- tids[[m]]
        train <- true_ft[tid, ]
        test <- true_ft[-tid, ]
        
        # fit parameters
        tabug_fit <- bn.fit(tabug, data = train, method = 'bayes')
        
        # assess performance
        predicts[[m]] <- t(attr(predict(tabug_fit, test, node = 'label',
                                        method = 'bayes-lw', prob = TRUE), 'prob'))[,2]
        labels[[m]] <-  test$label
}


## make plots ##
pred <- prediction(predicts, labels)
perf_roc <- performance(pred, "tpr","fpr")
perf_prc <- performance(pred, 'rec', 'prec')


dir.create(paste(data_type, 'perf', sep = "_"))
png(filename = 'perf/prc.png')
plot(perf_prc,col="grey82",lty=3)
plot(perf_prc,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE, main = 'roc', xlim = c(0,1), ylim = c(0,1))
dev.off()

png(filename = 'perf/roc.png')
plot(perf_roc,col="grey82",lty=3)
plot(perf_roc,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE, main = 'roc', xlim = c(0,1), ylim = c(0,1))
dev.off()