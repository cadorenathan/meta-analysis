
# read DGE files output from edgeR

DGE1 <- read.csv("DGE1.csv")
DGE2 <- read.csv("DGE2.csv")
DGE3 <- read.csv("DGE3.csv")

# create a list of the DGEs
DGE <- list(DGE1 = DGE1,
            DGE2 = DGE2,
            DGE3 = DGE3)


BiocManager::install("MetaVolcanoR", eval = FALSE)
library(MetaVolcanoR)


# performe the meta-analysis
meta_degs_comb <- combining_mv(diffexp=DGE,
                               pcriteria='FDR', 
                               foldchangecol='logFC',
                               genenamecol='SYMBOL',
                               geneidcol=NULL,
                               metafc='Mean',
                               metathr=0.01, 
                               collaps=TRUE,
                               jobname="MetaVolcano_Fisher",
                               outputfolder=".",
                               draw='HTML')


input <- meta_degs_comb@input
result <- meta_degs_comb@metaresult

write.csv(input, "meta_input.csv")
write.csv(result, "meta_result.csv")

