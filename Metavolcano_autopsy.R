
setwd("C:/Users/cador/OneDrive/Doutorado/Meta SBGM/Meta_autopsy_SVA")

GSE150316 <- read.csv("GSE150316_SVA.csv")
GSE155241 <- read.csv("GSE155241_SVA.csv")
GSE183533 <- read.csv("GSE183533_SVA.csv")

DGE <- list(GSE150316 = GSE150316,
            GSE155241 = GSE155241,
            GSE183533 = GSE183533)


#BiocManager::install("MetaVolcanoR", eval = FALSE)
library(MetaVolcanoR)

#esse d? o plot bonito
votecount <- MetaVolcanoR::votecount_mv(diffexp=DGE,
                                          pcriteria="FDR",
                                          foldchangecol='logFC',
                                          genenamecol='SYMBOL',
                                          geneidcol=NULL,
                                          pvalue=0.05,
                                          foldchange=0, 
                                          metathr=0.01,
                                          collaps=TRUE,
                                          jobname="MetaVolcano_VoteCount", 
                                          outputfolder=".",
                                          draw='HTML')

#esse ? o nosso
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

#esse s? se todos fossem microarray
meta_degs_rem <- rem_mv(diffexp=DGE,
                        pcriteria="FDR",
                        foldchangecol='logFC', 
                        genenamecol='Symbol',
                        geneidcol=NULL,
                        collaps=TRUE,
                        llcol='CI.L',
                        rlcol='CI.R',
                        vcol=NULL, 
                        cvar=TRUE,
                        metathr=0.01,
                        jobname="MetaVolcano",
                        outputfolder=".", 
                        draw='HTML',
                        ncores=1)


result_micro <- meta_degs_rem@metaresult
input_micro <- meta_degs_rem@input

#para ver um gene espec?fico
draw_forest(remres = meta_degs_rem,
             gene="SLC2A1",
             genecol="Symbol", 
             foldchangecol="logFC",
             llcol="CI.L", 
             rlcol="CI.R",
             jobname="MetaVolcano_plot",
             outputfolder=".",
             draw="PDF")

write.csv(result_micro, "result_micro.csv")
write.csv(input_micro, "input_micro.csv")

genes.map <-select(org.Hs.eg.db, as.character(result_micro$Symbol),c("SYMBOL","ENTREZID"), "SYMBOL")
colnames(genes.map)[1] <- "Symbol"
result_micro <- merge(result_micro, genes.map)

geneList <- result_micro$randomSummary
names(geneList) <- result_micro$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

gsea.go <- gseGO(geneList=geneList, 
                 ont ="ALL", 
                 keyType = "ENTREZID", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Hs.eg.db", 
                 pAdjustMethod = "fdr")
gsea.go.df <- gsea.go@result
gsea.go.sig <- gsea.go.df[which(gsea.go.df$p.adjust<0.05),]
gseaplot(gsea.go, by = "all", title = gsea.go$Description[1], geneSetID = 1)

write.csv(gsea.go.sig, "olf_gseaGO.csv")

gsea.kegg <- gseKEGG(geneList=geneList, 
                     organism = 'hsa',
                     nPerm = 10000, 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 1, 
                     verbose = TRUE, 
                     pAdjustMethod = "fdr")
gsea.kegg.df <- gsea.kegg@result
gsea.kegg.sig <- gsea.kegg.df[which(gsea.kegg.df$p.adjust<0.05),]
gseaplot(gsea.kegg, by = "all", title = gsea.kegg$Description[1], geneSetID = 1)


library("pathview")
hsa05171 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa05171",
                     species    = "hsa")

write.csv(gsea.kegg.sig, "olf_gseaKEGG.csv")
