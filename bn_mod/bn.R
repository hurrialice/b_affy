
## ----- misc -----------
setwd('~/b_affy/bn_mod')
# setwd('F:\\Project_group\\Git_Project\\FHIR\\b_affy\\bn_mod')

library(readr)
library(dplyr)
library(bnlearn)
library(Rgraphviz)
# shortcut function and load data
look <- function(m) m[1:10, 1:10]
# df containing filename, samplename and labels for each sample
ph <- read_rds('phenodata.rds')
# a table of differentially expressed genes
tab <- read_rds('degs.rds')
# data.matrix for each probe in each sample
data.matrix <- read_rds('data_matrix.rds')

## ------- notes ------------
cat('只是让脚本比较像R而且用kmeans做的binary离散化，后面建network的地方估计还要有一些tweak.大腿求抱emmm！')

## ------- functions --------------
# make raw features
make_feature_table <- function(ph. = ph, deg = tab, dm = data.matrix){
  sel.m <- dm[deg$PROBEID,]
  colnames(sel.m) <- ph.$sample_name[match(colnames(sel.m),ph.$file_name)]
  rownames(sel.m) <- deg$SYMBOL[match(rownames(sel.m),deg$PROBEID)]
  outdf <- t(sel.m) %>% as.data.frame()
  outdf$label <- ifelse(ph.$disease_state[match(rownames(outdf), ph.$sample_name)] == 'AMI', 1, 0)
  outdf
}


# make discrete
feature_cut <- function(ft = feature_table, nc = 2){
  ftm <- ft[,-which(names(ft) == "label")]
  #browser()
  #ftl <- split(ftm, col(ftm))
  bin_ft <- apply(ftm,2, function(f) {
    cldf <- kmeans(f, centers = nc, iter.max = 500, nstart = 10)
    ifelse(cldf$cluster == '1', 1, 0)
  })  %>% as.data.frame()
  #browser()
  bin_ft$label <- ft$label
  bin_ft
}


# featrue_filter 
featrue_filter <- function(df = discre_features, cut = 0.4){
  outcome <- df$label
  fcor <- apply(df, 2, function(f) cor(f, outcome, method = 'pearson'))
  sel.fs <- names(fcor)[abs(fcor) > cut]
  df[,sel.fs]
}

# partial integration
partial_integration <- function(dataset){
    # divide features into two/more groups
    # build two/more bayes networks with outcome
    # union of those networks
    return(union_network)
}


# structure learning
struc_learn <- function(dataset, algorithm='mmhc'){
  if(algorithm == 'mmhc'){
        net_stru = mmhc(dataset)
    } else if(algorithm == 'rsmax2'){
        net_stru = rsmax2(dataset)
    } else if(algorithm == 'tabu'){
        net_stru = tabu(dataset)
    } else if(algorithm == 'hc'){
        net_stru = hc(dataset)
    } else{
        stop(("algorithm name error !!"))
    }
  return(net_stru)
}

net_info <- function(net_stru){
  print(net_stru)
  bnlearn::graphviz.plot(net_stru, highlight = list(nodes = 'label'))
}

# parameter learning 
param_learn <- function(net_stru, dataset){
  cpt <- bn.fit(net_stru, dataset)
  return(cpt)
}


### ----- operations ---------
# get features

featrue_table <- make_feature_table()
discre_features <- feature_cut()
fin.fs <- featrue_filter()

# network learning

net_stru <- struc_learn(fin.fs, algorithm='mmhc')
cpts <- param_learn(net_stru, fin.fs)
net_info(net_stru)
