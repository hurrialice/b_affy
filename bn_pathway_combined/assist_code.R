## ------ data partitions ------
make_par <- function(){
  AMI_set <- make_ran(c(rep(1:10,2), 11))
  control_set <- make_ran(c(rep(1:10,3), 11))
  subjects <- c(AMI_set, control_set)
  caret::groupKFold(subjects, 11)
}

make_ran <- function(v) sample(v, length(v))

## ------- make pathway id2name conversion -------
make_path_dict <- function(hr = humanReactome){
  hr_id <- lapply(hr, function(p) p@id)
  names(hr_id) <- names(hr)
  return(hr_id)
}

getpathnames <- function(fc_tab){
  re_ids <- y$ID
  re_names <- names(pdict)[match(re_ids, pdict)]
}

## ------ discretilize ---------
make_discre <- function(m, psigma){
  o <- apply(m, 1, function(v){
    v0 <- (v - mean(v))/sd(v)
    out <- v0
    out[v0 < -psigma] <- -1
    out[abs(v0) < psigma] <- 0
    out[v0 > psigma] <- 1
    return(out)
  })
  t(o)
}

## ------- DE -----------
DE_pipe <- function(m, lab){
  # make matrix
  design <- model.matrix(~0 + lab)
  colnames(design) <- c("AMI", "control")
  fit <- lmFit(m, design)
  contrast.matrix <- makeContrasts(AMI-control, levels=design)
  
  # compare
  huvec_fits <- contrasts.fit(fit, contrast.matrix)
  huvec_ebFit <- eBayes(huvec_fits)
  tab <- topTable(huvec_ebFit, coef=1,lfc = 0.26, number = 1000)   
  tab <- tab[tab$P.Value < 0.05,]
  return(tab)
}


## ---------- entrez and symbols ---------
anno <- function(chip = hgu133plus2.db, t ){
  gs2 <- AnnotationDbi::select(x = chip, keys = rownames(tab), keytype = 'PROBEID',
                               columns = c('ENTREZID','GENENAME','SYMBOL')) %>% tbl_df()
  tab$PROBEID <- rownames(tab)
  tab <- left_join(tab, gs2) %>% tbl_df()
}


## --------- spia prepare -----------
spia_filter <- function(tab, id) {
  de_genes <- tab$ENTREZID
  all_path_nodes <- lapply(humanReactome, nodes)
  inter_path_nodes <- lapply(all_path_nodes, 
                             function(v) intersect(v, de_genes))
  
  all_path_nodes[lengths(inter_path_nodes) < 1] <- NULL
  interest_pathname <- names(all_path_nodes)
  path_subset <- humanReactome[interest_pathname]
  
  spia_rdata <- paste0('react_paths',id)
  
  prepareSPIA(path_subset, spia_rdata)
  de <- tab$logFC
  names(de) <- tab$ENTREZID
  all <- bg_genes$ENTREZID
  res <- runSPIA(de, all, spia_rdata)
  
  return(res)
}

## --------- FC prepare ------------
fc_filter <- function(tab){
  gl <- tab2genelist(tab)
  y <- gsePathway(gl, nPerm=1000,
                  minGSSize=10, pvalueCutoff=0.2,
                  pAdjustMethod="none", verbose=FALSE)
  y %>% tbl_df()
}

tab2genelist <- function(tab){
  gl <- tab$logFC
  names(gl) <- tab$ENTREZID
  sort(gl, decreasing = T)
}




## -------- make whitelist -------
# parse KGML file to graphNEL 
getg <- function(pathid){
  pathwayGraph(pathid)
  #url <- getKGMLurl(pathid, organism = "hsa")
  #g <- parseKGML2Graph(url)
}

# change names to entrez id

# entrez to symbol
entre2syms <- function(ev){
  tab$SYMBOL[match(ev, tab$ENTREZID)]
}
syms2entre <- function(ev){
  tab$ENTREZID[match(ev, tab$SYMBOL)]
}


# get connection between the connected nodes
# makeg <- function(adjm, spec_nodes, maxstep ){
#   g <- graph.adjacency(adjm)
#   nbl <- igraph::neighborhood(g, 
#                               order = maxstep, nodes = spec_nodes, mode = "out", mindist = 1)
#   nblist <- lapply(nbl, function(vs){
#     ids <- as_ids(vs)
#     intersect(spec_nodes, ids)
#   })
#   names(nblist) <- spec_nodes
#   dfout <-data.frame(from = rep(names(nblist), lengths(nblist)), to = unlist(nblist))
#   graph.data.frame(dfout)
# }

getsigpaths <- function(t){
  a <- t$Description[t$Description %in% names(pdict)]
  a[!is.na(a)]
}


make_wl <- function(t, a){
  
  sig.path.names <- getsigpaths(t)
  # print(sig.path.names)
  sig.paths <- humanReactome[sig.path.names] # as pathway list
  
  allglist <- lapply(sig.paths, getg)
  
  # browser()
  print('finished getg')
  # browser()
  allg_merg <- mergeGraphs(allglist)
  
  nonDE_genes <- bg_genes$ENTREZID[bg_genes$ENTREZID %in% tab$ENTREZID]
  # browser()
  merg_m <- as(allg_merg, "matrix")
  sel.geneids <- intersect(tab$ENTREZID, colnames(merg_m))
  nonsel.gids <- intersect(nonDE_genes, colnames(merg_m))
  misc <- colnames(merg_m)[!colnames(merg_m) %in% c(sel.geneids, nonsel.gids)]
  
  # better to make distinction
  align_vec <- c(sel.geneids, nonsel.gids, misc)
  merg <- merg_m[align_vec, align_vec]
  
  merg[sel.geneids, nonsel.gids] <- merg[sel.geneids, nonsel.gids]*a
  merg[nonsel.gids, sel.geneids] <- merg[nonsel.gids, sel.geneids]*a
  
  g <- graph.adjacency(merg, weighted = TRUE)
  
  # calculate distance of shortest paths
  dist_ioDE <- distances(g, v = sel.geneids, to = nonsel.gids, algorithm = "dijkstra")
  dist_ioDE2 <- distances(g, v = nonsel.gids, to = sel.geneids, algorithm = "dijkstra")
  dist_DE <- distances(g, v = sel.geneids, to = sel.geneids, algorithm = "dijkstra")
  
  dist_DE0 <- dist_DE[(!is.infinite(dist_DE))& dist_DE > 0]
  dist_ioDE[(!is.infinite(dist_ioDE))] -> dist_ioDE
  dist_ioDE2[(!is.infinite(dist_ioDE2))] -> dist_ioDE2
  dist2 <- c(dist_ioDE, dist_ioDE2)
  t.test(dist2, dist_DE0)
  
  dist_DE[dist_DE < a] <- 1
  dist_DE[dist_DE >= a] <- 0
  diag(dist_DE) <- 0
  deg <- graph.adjacency(dist_DE)
  

  # browser()
  # bet <- betweenness(deg)
  # degr <- degree(deg)
  return(deg)
  
  # # browser()
  # tonodes <- entre2syms(names(bet))[bet > 0 ] #& seq_along(bet) < 8]
  # df <- data.frame(from = 'label', to = tonodes)
  #return(degr)
}

## ----- get DE-based feature table ----
# mod_de_table <- function(){
#   #browser()
#   true_ft <- t(ds_ft[tab$PROBEID,tid]) %>% as.data.frame()
#   #browser()
#   colnames(true_ft) <- tab$SYMBOL[match(colnames(true_ft), tab$PROBEID)]
#   #browser()
#   true_ft <- cbind(true_ft, train_lab)
#   colnames(true_ft)[ncol(true_ft)] <- 'label'
#   print('fi mod')
#   o <- apply(true_ft, c(1,2), as.factor)
#   o <- as.data.frame(o)
#   
# }


## ----- use string protein interaction network ----

prepare_string <- function(tabs){
  dir.create("DE_genes")
  syms <- lapply(tabs, function(t) t$SYMBOL)
  for (i in seq_along(syms)){
    df <- data.frame(syms[[i]])
    write_csv(df, paste0('DE_genes/', i,'.csv'))
  }
}

getg_string <- function(path = 'stringnet_folds'){
  folds_name <- list.files(path)
  gs <- lapply(folds_name, function(nfold){
    tsv <- read_tsv(paste(path, nfold, sep = '/'))
    df <- data.frame(A = tsv$`#node1`, B = tsv$node2)
    graph.data.frame(df, directed = F)
  })
}
### ----- featrue tables --------
shrink_features <- function(ori_dm, fea_names, dict = bg_genes){
  m <- t(ori_dm)
  colnames(m) <- dict$SYMBOL[match(colnames(m), dict$PROBEID)]
  
  m_sub <- m[,fea_names]
  
  df <- as.data.frame(m_sub)
  df$label <- ifelse(ph$disease_state == 'AMI', 1, 0)
  df <- apply(df, c(1,2), as.factor)
  as.data.frame(df)
}

## ----- boot tabu -----------

# get the from-to dataframe
boot_tabu <- function(ori_tab){
  start_tan <- tree.bayes(ori_tab, 'label', head(colnames(ori_tab), -1))
  bt <- boot.strength(data = ori_tab, R = 500,algorithm = 'tabu', 
                      algorithm.args = list(start = start_tan, tabu = 50))
  #g <- averaged.network(bt, threshold = 0.7)
}


## ----- exact graphs --------
make_avgnet <- function(boot){
  bstr <- boot$strength
  bstr <- bstr[bstr > 0]
  bstr_cut_label <- max(boot$strength[boot$from == 'label' | boot$to == 'label']) - 0.00001
  bstr_cut_all <- quantile(bstr, 0.7)
  bstr_cut <- min(bstr_cut_all, bstr_cut_label)
  averaged.network(boot, threshold = bstr_cut_label)
}


mod_syms <- function(all, cs){
  two_comm <- unlist(cs[lengths(cs) <= 2])
  one_comm <- all[!all %in% unlist(cs)]
  
  #browser()
  cs[lengths(cs) <= 2] <- NULL
  ncom <- length(cs)
  names(cs) <- seq(ncom)
  cs[[ncom+1]] <- unname(c(one_comm, two_comm))
  names(cs) <- seq(ncom+1)
  return(cs)
}


getsubgraph <- function(arcdf){
  g <- graph.data.frame(arcdf)
  subcom <- subcomponent(g, 'label', 'all')
  sub_graph <- subgraph(g,subcom)
  as_edgelist(sub_graph)
}
