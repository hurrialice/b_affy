setwd('~/b_affy2/bn_mod')

library(readr)
library(dplyr)
library(igraph)
library(bnlearn)
library(Rgraphviz)
# original continous feature table
feature_table <- read_rds('featrue_table.rds')

## ----- global objective ------- 
# see how psigma permutation and feature_selection permutation can be realized in a possibily cleaner way.



## ----- sub-functions ----------
# split features to defined size
split_feature <- function(ft = feature_table, fea_groupSize = 15){
    label <- ft[,'label']
    ft <- ft[, -ncol(ft)]
    nfeatures <- ncol(ft)
    nflast <- nfeatures %% fea_groupSize
    
    ngroup <- ((nfeatures - nflast) / fea_groupSize) + 1
    group_id <- c(rep(seq(ngroup - 1), times = fea_groupSize), rep(ngroup, times = nflast))
    
    # lousy step, will simplify later
    sel_fea_tab <- lapply(lapply(seq(ngroup), function(i) which(group_id == i)), function(feaids) ft[,feaids])
    lapply(sel_fea_tab, function(m) cbind(m, label))
}



# permutate psigma, make it R-like
discrezation <- function(dfl, psigma){
  lapply(dfl, function(dataset){
    lb <- dataset[, which(colnames(dataset) == 'label')]
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
    #browser()
    out <- cbind(out, lb)
    out
  }
  )
}



struc_learn <- function(dataset, algorithm){
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

param_learn <- function(net_stru, dataset){
  cpt <- bn.fit(net_stru, dataset)
  return(cpt)
}

## -------- main-functions ----------
# permuate 2 params: featrue groups and psigma
make_pa2 <- function(ft, fea_groupSize){
  spfea <- split_feature()
  psigma_trials <- seq(from = 1, to = 10, by = 1)/10 # a small set (friendly for pc)
  mout <- mapply(discrezation, psigma = psigma_trials, MoreArgs = list(dfl = spfea))
  colnames(mout) <- psigma_trials
  mout
}



# select nodes and corresponding arcs that connects with 'lb'
get_arcs <- function(m, algo){
  mnets <- apply(m, c(1,2), function(l){
    dataset <- as.data.frame(l[[1]])
    
    # for each feature table, learn network
    net_stru <- struc_learn(dataset = dataset, algorithm = algo)
    arcdf <- net_stru$arcs
    
    # build igraph object
    g <- igraph::graph.data.frame(arcdf)
    
    if ('lb' %in% as_ids(V(g))){
      # order should be sufficiently large
      nbs <- neighborhood(g, order = 100, mode = 'all', mindist = 0, nodes = 'lb')
      connected_nodes <- as_ids(nbs[[1]])
      #browser()
      arcdf <- arcdf[arcdf[,1] %in% connected_nodes,]
  }
    else {arcdf <- matrix()}
    return(arcdf)}
    )
}



chose_one_psigma <- function(m, algo = 'max_nnodes', ftm){
  nnodes <- apply(m, c(1,2), function(l) length(l[[1]]))
  
  ### this chunk define how to chose psigma
  if (algo == 'max_nnodes'){
    max_nnodes<- apply(nnodes, 1, function(vec) unname(which(vec == max(vec)))[[1]])
  }
  else stop('algo not specified')
  ### however this method fails- -
  
  # get those subnets we chose
  subnets <- list()
  subfts <- list()
  for (i in seq(nrow(m))){
    j <- max_nnodes[i]
    subnets[[i]] <- m[i,j][[1]]
    subfts[[i]] <- ftm[i,j][[1]]
  }
  # remove empty net
  rm_ids <- which(sapply(subnets, function(m) length(m) == 1))
  subnets[rm_ids] <- NULL
  subfts[rm_ids] <- NULL
  #browser()
  nets <- do.call(rbind, subnets)
  nets <- nets[!duplicated(nets),]
  
  # get the features with their corresponding psigma discretized values
  nodes <- unique(unlist(nets))
  for (ii in seq_along(subfts)){
    subft <- subfts[[ii]]
    subfts[[ii]] <- subft[,colnames(subft) %in% nodes]
  }
  fts <- do.call(cbind, subfts)
  out <- list(features = fts, net = nets)
}





learn_params <- function(arc, ft){
  ft <- as.data.frame(ft[,!duplicated(colnames(ft))])
  e = empty.graph(colnames(ft))
  arcs(e) <- arc
  param_learn(net_stru = e , dataset = ft)
}


## ------ operations -----------
pa2 <- make_pa2(ft = feature_table, fea_groupSize = 15)

related_arcs <- get_arcs(m = pa2, algo = 'mmhc')
#!! warnings from bnlearn, problem unidentified!... big leg please check here..QaQ

subarcs <- chose_one_psigma(m = related_arcs, algo = 'max_nnodes', ftm = pa2)

cpts <- learn_params(arc = subarcs$net, ft = subarcs$features)


