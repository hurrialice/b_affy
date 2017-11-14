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
