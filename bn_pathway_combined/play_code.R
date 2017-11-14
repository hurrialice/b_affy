chip <- hgu133plus2.db
bg_genes <- AnnotationDbi::select(x = chip, keys = rownames(data.matrix), keytype = 'PROBEID',
                                  columns = c('ENTREZID','GENENAME','SYMBOL')) %>% tbl_df()

DE_genes <- tab$ENTREZID
nonDE_genes <- bg_genes$ENTREZID[!bg_genes$ENTREZID %in% DE_genes]


t <- fc_res[[1]]
sig.path.names <- getsigpaths(t)
# print(sig.path.names)
sig.paths <- humanReactome[sig.path.names] # as pathway list

allglist <- lapply(sig.paths, getg)

# browser()
allg_merg <- mergeGraphs(allglist)
# browser()
merg_m <- as(allg_merg, "matrix")
sel.geneids <- intersect(tab$ENTREZID, colnames(merg_m))
nonsel.gids <- intersect(nonDE_genes, colnames(merg_m))
misc <- colnames(merg_m)[!colnames(merg_m) %in% c(sel.geneids, nonsel.gids)]

align_vec <- c(sel.geneids, nonsel.gids, misc)
merg <- merg_m[align_vec, align_vec]

a <- 5
merg[sel.geneids, nonsel.gids] <- merg[sel.geneids, nonsel.gids]*a
merg[nonsel.gids, sel.geneids] <- merg[nonsel.gids, sel.geneids]*a


g <- graph.adjacency(merg, weighted = TRUE)
dist_ioDE<- distances(g, v = sel.geneids, to = nonsel.gids, algorithm = "dijkstra")
dist_ioDE2 <- distances(g, v = nonsel.gids, to = sel.geneids, algorithm = "dijkstra")
dist_DE <- distances(g, v = sel.geneids, to = sel.geneids, algorithm = "dijkstra")

dist_DE0 <- dist_DE[(!is.infinite(dist_DE))& dist_DE > 0]
dist_ioDE[(!is.infinite(dist_ioDE))] -> dist_ioDE
dist_ioDE2[(!is.infinite(dist_ioDE2))] -> dist_ioDE2
dist2 <- c(dist_ioDE, dist_ioDE2)
t.test(dist2, dist_DE0)

dist_DE[dist_DE < 5] <- 1
dist_DE[dist_DE >= 5] <- 0
diag(dist_DE) <- 0
deg <- graph.adjacency(dist_DE)
##############

tabu_boot <- read_rds('tabugs_new.rds')
tabu_boot <- lapply(tabu_boot, function(t) t[t[,'direction'] > 0.5,])



#############

# use hybrid method here
el_gs <- list()
for (i in seq_along(string_gs)){
  tid <- tids[[i]]
  tab <- tab_results[[i]]
  g <- string_gs[[i]]
  comm <- cluster_edge_betweenness(g)
  c_syms <- communities(comm)
  
  css <- mod_syms(all = tab$SYMBOL, cs = c_syms)
  
  # init containers
  boots <- list()
  subgs <- list()
  rel.edges <- list()
  
  for (j in seq_along(css)){
    cs <- intersect(css[[j]], bg_genes$SYMBOL)
    ft <- shrink_features(ori_dm = ds_ft[tid,], fea_names = cs, dict = bg_genes)
    boot <- boot.strength(ft, R = 200, cpdag = FALSE, 
                               algorithm = 'mmhc', 
                               algorithm.args = list(maximize.args = 
                                                       list(start = naive.bayes(ft,'label'))))
    print('bootstrape done')
    boot <- boot[boot$direction > 0.5,]
    boots[[j]] <- boot
    
    subgs[[j]] <- make_avgnet(boot)
    
    arc_sets <- arcs(subgs[[j]])
    rel.edges[[j]] <- getsubgraph(arc_sets) 
  }
  el_gs[[i]] <- do.call(rbind, rel.edges)
}

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

