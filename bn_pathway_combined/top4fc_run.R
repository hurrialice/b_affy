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

source('otherfs_assist.R')

## this time we choose the top 4 features with the most fold change

# ----- prepare ------------
data.matrix <- read_rds('raw1106/data_matrix.rds')
ph <- read_rds('raw1106/phenodata.rds')
ds_ft <- make_discre(data.matrix, 0.6)
# map to chip
chip <- hgu133plus2.db
bg_genes <- AnnotationDbi::select(x = chip, keys = rownames(data.matrix), keytype = 'PROBEID',
                                  columns = c('ENTREZID','GENENAME','SYMBOL')) %>% tbl_df()
# ----- create data partitions --------
set.seed(1234)
tids <- make_par()

# ----- start with tab_list ----------- ### further change here!
tabs <- read_rds('tabs_1109.rds')
# enter manually
n <- 12

# ------ operations -------------------
dm <- make_discre(data.matrix, 0.6)
tabugs <- read_rds('tabugs.rds')
nfeatures <- lapply(tabugs, function(bn) length(mb(bn, 'label')))

# ------- functions -----------------
for (i in seq_along(tids)){
  tid <- tids[[i]]
  tab <- tabs[[i]]
  n <- nfeatures[[i]]
  # remain only the selected features
  features <- topNfs(n, tab)
  print(features)
  true_ft <- shrink_features(ori_dm = dm, fea_names = features)
  
  # contrast train from test
  train_ft <- true_ft[tid,]
  test_ft <- true_ft[-tid,]
  
  # build network
  tan <- tree.bayes(train_ft, 'label', head(colnames(true_ft), -1))
  tan.fit <- bn.fit(tan, data = true_ft, method = 'bayes')
  
  # test performance
  predicts[[i]] <- t(attr(predict(tan.fit, test_ft, prob = TRUE), 'prob'))[,2]
  labels[[i]] <-  test_ft$label
}

dir.create('top4_perf')
pred <- prediction(predicts, labels)
perf_roc <- performance(pred, "tpr","fpr")
perf_auc <- performance(pred, 'auc')
perf_prc <- performance(pred, 'rec', 'prec')

png(filename = 'top4_perf/prc_rtx.png')
plot(perf_prc,col="grey82",lty=3)
plot(perf_prc,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE, main = 'roc', xlim = c(0,1), ylim = c(0,1))
dev.off()

png(filename = 'top4_perf/roc_rtx.png')
plot(perf_roc,col="grey82",lty=3)
plot(perf_roc,lwd=3,avg="vertical",spread.estimate="boxplot",add=TRUE, main = 'roc', xlim = c(0,1), ylim = c(0,1))
dev.off()

auc <- unlist(perf_auc@y.values) %>% mean
