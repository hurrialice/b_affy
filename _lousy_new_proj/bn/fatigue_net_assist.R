library(bnlearn)
library(readr)
library(dplyr)
library(caret)
library(limma)
library(ROCR)
library(future)
library(annotate)
library(hgu133plus2.db)

# get background genes
chip <- hgu133plus2.db
bg_genes <- read_rds("bg_genes.rds")
#bg_genes <- AnnotationDbi::select(x = chip, keys = rownames(data.matrix), keytype = 'PROBEID',
#                                 columns = c('ENTREZID','GENENAME','SYMBOL')) %>% tbl_df()


# for the sake of easier visualization
look <- function(m) m[1:10, 1:10]

# plan for parallel processing via future
future::plan(multicore)

## ------ data partitions ------
set.seed(1234)
make_par <- function(){
        nsamples <- ncol(data.matrix)
        trainid <- list()
        # browser()
        for (i in 1:nsamples){
                a <- as.integer(i)
                # note long vector be placed first
                trainid[[i]] <- setdiff(seq(nsamples), a)
        }
        return(trainid)
}



# init list-like object for data storage
init_list <- function(clist){
        for (i in clist){
                assign(x = i, value = list(), envir = globalenv())
        }
}

## ------ discretilize by gene ---------
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
        tab <- topTable(huvec_ebFit, coef=1,lfc = 0.4, number = 1000)   
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



### ----- featrue tables --------
shrink_features <- function(ori_dm, fea_names, dict = bg_genes){
        m <- t(ori_dm)
        m_sub <- m[,fea_names]
        #browser()
        colnames(m_sub) <- dict$SYMBOL[match(colnames(m_sub), dict$PROBEID)]
        
        df <- as.data.frame(m_sub)
        df$label <- ifelse(ph$disease_state == 'AMI', 1, 0)
        df <- apply(df, c(1,2), as.factor)
        as.data.frame(df)
}

## ----- boot tabu -----------

# get the from-to dataframe
boot_tabu <- function(ori_tab){
        start_tan <- tree.bayes(ori_tab, 'label', head(colnames(ori_tab), -1))
        bt <- boot.strength(data = ori_tab, R = 200, m = nrow(ori_tab), algorithm = 'tabu', debug = T,
                            algorithm.args = list(start = start_tan, tabu = 50))
}
