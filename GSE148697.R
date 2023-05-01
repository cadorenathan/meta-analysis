
#install package
if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("edgeR")
n

library("edgeR")

#select folder
setwd("C:/Users/cador/OneDrive/Doutorado/Meta SBGM/GSE148697")

#to check the folder
getwd()

#first thing is to call my reading count files
#they are separate files, so I must use the readDGE function
#this function will not only read, but collate them all into a single file
#if the files were already joint, I should use read.delim function

files <- c("GSM4476798.tabular",
           "GSM4476799.tabular",
           "GSM4476800.tabular",
           "GSM4476801.tabular",
           "GSM4476802.tabular",
           "GSM4476803.tabular")
DG <- readDGE(files,header=F)

#readDGE function appeared in the Global Environment
#then, DG value appeared with the 7 files I indicated

#need to create a data frame now

# load each sample individually
covid1 = read.table('GSM4476798.tabular',sep='\t',header=FALSE)
covid2 = read.table('GSM4476799.tabular',sep='\t',header=FALSE)
covid3 = read.table('GSM4476800.tabular',sep='\t',header=FALSE)

control1 = read.table('GSM4476801.tabular',sep='\t',header=FALSE)
control2 = read.table('GSM4476802.tabular',sep='\t',header=FALSE)
control3 = read.table('GSM4476803.tabular',sep='\t',header=FALSE)

# combine data sets into a matrix
geneCounts = data.frame(control1[,2], control2[,2],control3[,2],
                        covid1[,2], covid2[,2], covid3[,2])
row.names(geneCounts) = control1[,1]
sizeGeneCounts = dim(geneCounts)
geneCounts = geneCounts[1:(sizeGeneCounts[1]-5),]

#follow genecounts order of conditions
condition = c(rep('control',3),rep('covid',3))
sampleNames = c('control1','control2','control3',
                'covid1','covid2','covid3')
colnames(geneCounts) = sampleNames
View(geneCounts)

# build the generalized linear model that will be used for differential expression testing
dge <- DGEList(counts=geneCounts, group=condition)
design <- model.matrix(~condition+0, data=dge$samples)
colnames(design) = gsub("condition","",colnames(design))

# perform normalization by TMM
dge <- calcNormFactors(dge)
plotMDS(dge)
norm_counts <- cpm(dge,log = TRUE, prior.count = 3)
exp <- as.data.frame(norm_counts)

#write plot
jpeg(file="GSE148697_MDSplot.jpeg", width=5000, height=5000, units="px", res=300)
plotMDS(dge)
dev.off()

#estimate dispersion
disp <- estimateGLMCommonDisp(dge, design)
disp <- estimateGLMTrendedDisp(disp, design)
disp <- estimateGLMTagwiseDisp(disp, design)
plotBCV(disp)

#write plot
jpeg(file="GSE148697_dispersion_plot.jpeg", width=5000, height=5000, units="px", res=300)
plotBCV(disp)
dev.off()

########## PCA ANALYSIS #############

#taken from microarray assays

## load the exp file
geneCounts_PC<-geneCounts

## expA_PC must a be a matrix
geneCounts_PC<-data.matrix(geneCounts_PC)

## transpose the expPC
transp<-t(geneCounts_PC)

## analyse principal component
p3A<-prcomp(transp, retx=TRUE, center=TRUE, scale=FALSE)
summary(p3A)

## associate to p3scores the values of PC
p3Ascores <- p3A$x

######################################################################

##### PCA PLOT #####

#also taken from PCA analysis

jpeg(file="GSE148697_PCA.jpeg", width=3200, height=3200, units="px", res=300)

## plot the PC
plot(p3Ascores[,1], p3Ascores[,2], xlab="PCA 1", ylab="PCA 2",
     type="p", cex.lab=0.75, cex.axis=0.75, 
     #xlim=c(-200,250), ylim=c(-200,170),
     col=c('deepskyblue','deepskyblue','deepskyblue','darkred','darkred',
           'darkred'),
     main="PCA scores", cex.main=1.2, font.main=1,pch=15)

## apply to plot
text(p3Ascores, colnames(sampleNames), cex=0.5, pos=4, col="black")

legend("topleft", legend=c("Control","COVID"),
       bty="n", xjust = 1, yjust = 1,
       cex=.75, y.intersp=1, col=c('deepskyblue', 'darkred'), pch=20)

##
dev.off()

#######################################################################
##### DIFFERENTIAL EXPRESSION ANALYSIS #####

fit <- glmFit(disp, design)
lrt <- glmLRT(fit)
topTags(lrt)

#at last ¬¬

# tell edgeR what comparison you want to perform
covidVScontrol = makeContrasts(covid-control, levels=design)

# obtain gene symbol
BiocManager::install("org.Hs.eg.db")
n

library(org.Hs.eg.db)
genes.map <-select(org.Hs.eg.db, as.character(control1$V1),c("SYMBOL","ENTREZID"), "ENTREZID")

# perform the differential expression testing for that comparison
lrt.covidVScontrol = glmLRT(fit, contrast=covidVScontrol)
res.covidVScontrol<-topTags(lrt.covidVScontrol, n=60000, sort.by = "p.value")

#add gene names
table.covidVScontrol <- as.data.frame(res.covidVScontrol$table)
table.covidVScontrol$ENTREZID <- row.names(table.covidVScontrol)
table.covidVScontrol <- merge(table.covidVScontrol,genes.map)

#creating file
write.csv(table.covidVScontrol, file="GSE148697.csv")

#aqui no script p.value é ajustado
#PValue não ajustado


