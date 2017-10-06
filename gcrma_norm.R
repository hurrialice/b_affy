library(GEOquery)
library(affy)
library(dplyr)
library(readr)


## --- get files ---------
getGEOSuppFiles("GSE48060")
untar("GSE48060/GSE48060_RAW.tar", exdir="data")
cels <- list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep="/"), gunzip)
# set this path to which cotains the soft file, [gz] is acceptable
gds5074 <- getGEO(filename='GDS5074.soft.gz')




#### ---- make labels --------
# basic info: platform, 
Meta(gds5074)
# data
Table(gds5074)
# exact the labels from soft files
eset <- GDS2eSet(gds5074)
pheno_raw <- pData(eset) %>% tbl_df()
# make phenodata object
files <- system("ls data", intern = TRUE)
mapply(function(p,x) grepl(p, x),
       p = pheno_raw$sample,
       x = files) 
# problem identified, require careful curation because the naming sequence is not identical
ph <- tibble(sample_name = pheno_raw$sample,
             file_name = files[mapply(function(p,x) grep(p,x),
                                      p = pheno_raw$sample,MoreArgs = list(x = files))],
             disease_state = c(rep('AMI', 31), rep("control", 21)),
             recur = NA) %>% dplyr::select(file_name, sample_name, disease_state, recur)

# modify the labels
ph$recur[grepl('without', pheno_raw$description)] <- 'AMI_no'
ph$recur[grepl('control', pheno_raw$description)] <- 'control'        
ph$recur[is.na(ph$recur)] <- 'AMI_yes'
# have a look
ph
write_tsv(ph, 'data/phenodata.txt')




## ----- gcrma normalize ----------
library  (simpleaffy)
celfiles <- read.affy(covdesc="phenodata.txt", path="data")     # affybatch object
celfiles.gcrma <- gcrma(celfiles)                               # expressionset object for normalized data
celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=TRUE, remove.dupEntrez=FALSE)
celfiles.filtered$filter.log

data.matrix <- celfiles.filtered$eset %>% exprs





## ---- clustering ---------
# dimmension reduction using PCA
color=c(rep('red', 31), rep('green', 21))
data.PC = prcomp(t(data.matrix),scale.=FALSE)
PCdf <- data.PC$x %>%  as.data.frame() 
pairs(~PC1+PC2+PC3,data=PCdf, 
      main="PCA", col = color)
PCdf$label <- c(rep('AMI', 31), rep('control', 21)) %>% as.factor


# dimmension reduction using tsne
set.seed(1)
data.tsne <-  Rtsne(t(data.matrix), dims = 2, perplexity=10, verbose=TRUE, max_iter = 500)

colors = rainbow(length(unique(train$label)))
names(colors) = unique(train$label)
plot(data.tsne$Y, t='n', main="tsne")
text(data.tsne$Y, labels=train$label, col=colors[train$label])





# make it interactive using plotly
library(plotly)
Sys.setenv("plotly_username"="hurrialice")
Sys.setenv("plotly_api_key"="CYOeIoOFiGxEFy8sgcwc")


p <- plot_ly(PCdf, x = ~PC1, y = ~PC2, z = ~PC3, color = ~label, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))
chart_link = api_create(p, filename="PCA_affy")
chart_link

set.seed(1)
data.tsne3 <-  Rtsne(t(data.matrix), dims = 3, perplexity=10, verbose=TRUE, max_iter = 500)
tsne_df <- as.data.frame(data.tsne3$Y)
colnames(tsne_df) <- c('tsne_1', 'tsne_2', 'tsne_3')
tsne_df$label <- PCdf$label
p2 <- plot_ly(tsne_df, x = ~`tsne_1`, y = ~`tsne_2`, z = ~`tsne_3`, color = ~label, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'TSNE1'),
                      yaxis = list(title = 'TSNE2'),
                      zaxis = list(title = 'TSNE3')))
chart_link2 = api_create(p2, filename="tSNE_affy")
chart_link2


# hclust
library(dendextend)
distance <- dist(PCdf[,-1],method="maximum")
dendo <- hclust(distance) %>% as.dendrogram()
new_label <- ph$sample_name[match(labels(dendo), ph$file_name)]
new_phl <- ph$recur[match(labels(dendo), ph$file_name)]
new_phl <- gsub('control','green', new_phl)
new_phl[!grepl('green', new_phl)] <- 'red'


dendo %>% set('labels', new_label) %>% 
  set("labels_col", new_phl) %>% # change color
  set("labels_cex", 1) %>% # Change size
  plot(main = "hierarchical clustering of samples") # plot




## ------ limma ------------
# make design matrix
samples <- celfiles.gcrma$disease_state %>% as.factor()
design <- model.matrix(~0 + samples)
colnames(design) <- c("AMI", "control")
# use limma
library(limma)
fit <- lmFit(exprs(celfiles.filtered$eset), design)
contrast.matrix <- makeContrasts(AMI-control, levels=design)
huvec_fits <- contrasts.fit(fit, contrast.matrix)
huvec_ebFit <- eBayes(huvec_fits)
tab <- topTable(huvec_ebFit, coef=1,lfc = 0.5, number = 1000)   



  
## ----- annotate by gene symbols ---------
library(hgu133plus2.db)
library(annotate)
chip <- hgu133plus2.db
gs2 <- AnnotationDbi::select(x = chip, keys = rownames(tab), keytype = 'PROBEID',
                             columns = c('ENTREZID','GENENAME','SYMBOL')) %>% tbl_df()
tab$PROBEID <- rownames(tab)
tab <- left_join(tab, gs2) %>% tbl_df()
write_csv(tab, 'enrich_tab.csv')



# ----- reactome api -------------
require(ReactomePA)
x <- enrichPathway(gene=tab$ENTREZID, pvalueCutoff=0.05, readable=T, pAdjustMethod = 'none') %>% as.data.frame()
# abort! maybe api too old? not enough meaningful result..
# use reactome-web derived csv
react <- read_csv('reactome_out.csv') %>% dplyr::filter(`Entities pValue` < 0.05)
react_shunk <- react %>% dplyr::select(pathway = `Pathway name`, pVal = `Entities pValue`)
write_csv(react_shunk, 'reactome_sig.csv')



## ----- gprofiler ------------
# get enriched GO terms and KEGG terms
library(gProfileR)
gdf <- gprofiler(tab$SYMBOL) %>% tbl_df()
write_csv(gdf, 'gprofiler_enriched.csv')


## ----- bicluster heatmap -------
library(d3heatmap)
tab0 <- tab[ order(-tab$logFC), ] %>% dplyr::filter(P.Value < 0.0005)
de_pids <- tab0$PROBEID
dem <- data.matrix[de_pids,]
colnames(dem) <- ph$sample_name[match(colnames(dem), ph$file_name)]
rownames(dem) <- tab$ENTREZID[match(rownames(dem), tab$PROBEID)]
d3heatmap(dem, scale = "col", dendrogram = "row" )
