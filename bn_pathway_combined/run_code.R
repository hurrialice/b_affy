#setwd("/home/qingzhang/b_affy/bn_mod")

# load library
library(bnlearn)
library(caret)
library(graphite)
library(readr)
library(dplyr)
library(igraph)
library(KEGGgraph)
library(KEGG.db)
library(SPIA)
library(limma)
library(hgu133plus2.db)
library(annotate)
library(ROCR)
library(future)
library(ReactomePA)
plan(multicore)
look <- function(m) m[1:5, 1:5]
# get functions
source('assist_code.R')
# load data set
data.matrix <- read_rds('raw1106/data_matrix.rds')
ph <- read_rds('raw1106/phenodata.rds')
dir_name <- 'gsea_paths'
humanReactome <- read_rds('humanReactome.rds')
pdict <- make_path_dict()

# make discrete data matrix for all genes
ds_ft <- make_discre(data.matrix, 0.6)

# create the container for output in this run
dir.create(dir_n)


# # create the identical partition of sample
# set.seed(1234)
# tids <- make_par()
# write_rds(tids, 'tids_uniform.rds')


# get background gene
chip <- hgu133plus2.db
bg_genes <- AnnotationDbi::select(x = chip, keys = rownames(data.matrix), keytype = 'PROBEID',
                             columns = c('ENTREZID','GENENAME','SYMBOL')) %>% tbl_df()

## init containers
# tab_results <- list()
spia_ids <- list()
path_graph_ids <- list()
tabugs <- list()
tabugsf <- list()

predicts <- list()
labels <- list()
# spia_resf <- list()
fc_resf <- list()
true_fts <- list()

train_labels <- list()1
2
3
4
data(learning.test)
res = gs(learning.test)
cpdag(res)
vstructs(res)
test_labels <- list()
# ## ---- load resumed files -----
# tids <- read_rds('tids_uniform.rds')
# spia_res <- read_rds('spia_tables.rds')
# tab_results <- read_rds('tabs_1109.rds')
## the real exciting cycle begins = =

# ---------- DE --------------------
for (i in seq_along(tids)){

  # define train and test
  tid <- tids[[i]]
  train <- data.matrix[,tid]
  test <- data.matrix[, -tid]

  # differential expression
  train_label <- as.factor(ph$disease_state[tid])
  tab <- DE_pipe(train, train_label)

  # give labels to DE genes
  tab <- anno(t = tab)
  tab_results[[i]] <- tab
  print('tab result saved')
  print(tab)

  # modify label
  train_labels[[i]] <- ifelse(train_label == 'AMI', 1, 0 )
  test_labels[[i]] <- ifelse(as.factor(ph$disease_state[-tid]) == 'AMI', 1, 0)
  
  # real feature table for tabu search (train + test)
  true_fts[[i]] <- shrink_features(ori_dm = ds_ft, fea_names = tab$SYMBOL)
}

# save: tabs
write_rds(tab_results, 'tabs_1114.rds')
warnings()

# ## -------- pathways enrichment --------------------
# for (j in seq_along(tab_results)){
#   tab <- tab_results[[j]]
#   # spia_resf[[j]] <- future({spia_filter(tab, j)})
#   fc_resf[[j]] <- future({fc_filter(tab)})
# }
# 
# 
# fc_res <- lapply(fc_resf, FUN = value)
# write_rds(fc_res, 'fc_res.rds')
# print('FC finished!!')
# 
# # see pathways
# 
# fc_overlap  <- unlist(lapply(fc_res, function(t) t$Description)) %>% table %>% sort()
# spia_overlap <- unlist(lapply(spia_res, function(t) t$Name[t$pGFdr < 0.1])) %>% table %>% sort()
# 
# 
# # save SPIA
# # write_rds(spia_res, 'spia_tables.rds')
# 
# # for (m in seq_along(spia_res)){
# #   spia_tab <- spia_res[[m]]
# #   sig.pathids_1 <- spia_tab$Name[spia_tab$pGFdr < 0.1 | spia_tab$pGFWER < 0.1]
# #   sig.pathids_2 <- spia_tab$Name[1]
# #   sig.pathids <- unique(c(sig.pathids_2, sig.pathids_1)) # at least one pathway will remain
# #   spia_ids[[m]] <- sig.pathids
# # }
# 
# bets <- list()
# 
# for (k in seq_along(fc_res)){
#   fc_tab <- fc_res[[k]]
#   tab <- tab_results[[k]]
#   bets[[k]] <- make_wl(fc_tab, a = 5)
#   #path_graph_ids[[k]] <- wl$to
# }
# # save whitelist
# write_rds(path_graph_ids, 'wl_tonodes.rds')
# wls <- lapply(path_graph_ids, function(d){
#   data.frame(from = 'label', to = d)
# })
# print('whitelist made!')
# write_rds(true_fts, 'true_fts.rds')
# print('fts made')


# input prior knowledge from STRING database
prepare_string(tab_results)
string_gs <- getg_string()


## ------ make tabu -------------
for (mm in seq_along(true_fts)){
  true_ft <- true_fts[[mm]]
  tid <- tids[[mm]]
  train <- true_ft[tid, ]
  wl <- wls[[mm]]
  # boot tabu, return a graph
  tabugsf[[mm]] <- future({boot_tabu(wl, ori_tab = train)})
}

tabugs <- lapply(tabugsf, FUN = value)
write_rds(tabugs, 'tabugs.rds')
paste0('tabu finished at ', Sys.time())


tabugs <- read_rds('tabugs.rds')
true_fts <- read_rds('true_fts.rds')


## ------- final assess -----------
for (m in seq_along(tabugs)){
  tabug <- tabugs[[m]]
  true_ft <- true_fts[[m]]
  tid <- tids[[m]]
  train <- true_ft[tid, ]
  test <- true_ft[-tid, ]
  
  # make final TAN
  sig.nodes <- mb(tabug, 'label')
  fin_ft <- train[,c(sig.nodes, 'label')]
  final_tan <- tree.bayes(x = fin_ft, 'label', sig.nodes)
  tan.fit <- bn.fit(final_tan, data = fin_ft, method = 'bayes')
  
  # assess performance
  predicts[[m]] <- t(attr(predict(tan.fit, test, prob = TRUE), 'prob'))[,2]
  labels[[m]] <-  test$label
}
save(tids, tab_results, spia_ids, path_graph_ids, tabugs, predicts, file = 'exact_future.RData')

pred <- prediction(predicts, labels)
perf_roc <- performance(pred, "tpr","fpr")
perf_auc <- performance(pred, 'auc')
perf_prc <- performance(pred, 'rec', 'prec')


dir.create('perf')
png(filename = 'perf/prc_r1.png')
plot(perf_prc,col="grey82",lty=3)
plot(perf_prc,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE, main = 'roc', xlim = c(0,1), ylim = c(0,1))
dev.off()

png(filename = 'perf/roc_r1.png')
plot(perf_roc,col="grey82",lty=3)
plot(perf_roc,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE, main = 'roc', xlim = c(0,1), ylim = c(0,1))
dev.off()

save(tab_results, spia_ids, path_graph_ids, tabugs, file = 'exact1107.RData')


sel.features <- lapply(tabugs, function(bn) mb(bn, 'label'))
