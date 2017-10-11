look <- function(m) m[1:10, 1:10]
ph <- read_tsv('data/phenodata.txt')


library(readr)
selected.probes <- tab$PROBEID
probe_exp <- data.matrix[selected.probes,]
rownames(probe_exp) <- tab$SYMBOL
bndf0 <- t(probe_exp)
rownames(bndf0) <- ph$sample_name[match(rownames(bndf0), ph$file_name)]
bndf <- as.data.frame(bndf0) #%>% dplyr::select(-one_of('NASP','GOLGA8N','ADAMTS1','PRMT6','PPM1L'))

bndf <- discretize(bndf,method = 'quantile', breaks = 3)
