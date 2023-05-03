# BiocManager::install("DOSE")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("KEGGREST")
# BiocManager::install("ReactomePA")

# Load packages for functional enrichment

library("DOSE")
library("clusterProfiler")
library("org.Hs.eg.db") # for ids mapping
library("KEGGREST")

# load your gene list and atribute to "genes.map" 

## perform enrichment for all target genes found (using GO)
ego <- enrichGO(gene          = unique(na.omit(genes.map$ENTREZID)),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05)
enrich.go<- setReadable(ego,org.Hs.eg.db,keyType = "ENTREZID")
go.table <- as.data.frame(enrich.go)

#drop GO level 4
ego2<-dropGO(ego,level=6)

ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)
go2.table <- as.data.frame(ego2)
View(go2.table)

write.csv(go.table, "GO_targets.csv")
write.csv(go2.table, "GO2_table.csv")

jpeg(file="dotplot_KEGG.jpeg", width=6200, height=4800, units="px", res=600)
dotplot(enrich.kegg, showCategory=20, orderBy="GeneRatio")
dev.off()

jpeg(file="dotplot_targets.jpeg", width=6200, height=4800, units="px", res=600)
dotplot(ego, showCategory=20)
dev.off()
enrichMap(ego)

jpeg(file="cnetplot_GO.jpeg", width=6200, height=4800, units="px", res=600)
cnetplot(ego, categorySize="pvalue", showCategory = 3, fixed = TRUE, node_label='all')
dev.off()

##############################################################################

library("ReactomePA")

rPA <- enrichPathway(gene=genes.map$ENTREZID,pvalueCutoff=0.1, readable=T)
rPA2 <- as.data.frame(rPA)
View(rPA2)
write.csv(rPA2, "Reactome_results.csv")

jpeg(file="dotplot_Reactome.jpeg", width=6200, height=4800, units="px", res=600)
dotplot(rPA, showCategory=20)
dev.off()


#####################################################

genes.code <- genes.map$ENTREZID

#enrich disease ontology
enrich.dis <- enrichDO(gene          = genes.code,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(genes.code),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

disease.read <- setReadable(enrich.dis,'org.Hs.eg.db')
disease.table <- disease.read@result


#enrich disGenNet
dgn <- enrichDGN(genes.map$ENTREZID)
dgn.read <- setReadable(dgn,'org.Hs.eg.db')
dgn.table <- dgn.read@result

#plots
library(enrichplot)

jpeg(file="dotplot_DGN.jpeg", width=8200, height=4800, units="px", res=600)
dotplot(dgn, showCategory=20)
dev.off()

write.csv(dgn.table, "dgn_results.csv")

########################################################

# perform enrichment for all target genes found (using KEGG.db)
enrich.kegg <- enrichKEGG(gene = unique(na.omit(genes.map$ENTREZID)),
                          organism     = "human",
                          pvalueCutoff = 0.1,
                          use_internal_data = FALSE,
                          pAdjustMethod = "fdr")

enrich.kegg <- setReadable(enrich.kegg,org.Hs.eg.db,keyType = "ENTREZID")
kegg.table <- as.data.frame(enrich.kegg)
View(kegg.table)
write.csv(kegg.table, "KEGGresults.csv")

# browseKEGG(enrich.kegg, 'hsa04060')


