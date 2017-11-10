
topNfs <- function(n, tab){
  logfc <- tab$logFC
  names(logfc) <- tab$SYMBOL
  sorted <- sort(abs(logfc), decreasing = T)[1:n]
  names(sorted)
}

ranNf <- function(n, tab){
  sample(tab$SYMBOL, size = n)
}


# make true feature table
shrink_features <- function(ori_dm, fea_names, dict = bg_genes){
  m <- t(ori_dm)
  colnames(m) <- dict$SYMBOL[match(colnames(m), dict$PROBEID)]
  
  m_sub <- m[,fea_names]
  
  df <- as.data.frame(m_sub)
  df$label <- ifelse(ph$disease_state == 'AMI', 1, 0)
  df <- apply(df, c(1,2), as.factor)
  as.data.frame(df)
}


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
