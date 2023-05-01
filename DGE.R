
# install package
 if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("edgeR")

library("edgeR")

# to call the reading count files output from featureCounts tool (Galaxy Server)
# they are separate files, so you must use the readDGE function
# this function will not only read, but collate them all into a single file
# if the files were already joint, you should use read.delim function

files <- c("file1.tabular",
           "file2.tabular",
           "file3.tabular",
           "file4.tabular",
           "file5.tabular",
           "file6.tabular")
DG <- readDGE(files,header=F)

# readDGE function appeared in the Global Environment
# then, DG value appeared with the 6 files I indicated

# need to create a data frame now

# load each sample individually
covid1 = read.table('file4.tabular',sep='\t',header=FALSE)
covid2 = read.table('file5.tabular',sep='\t',header=FALSE)
covid3 = read.table('file6.tabular',sep='\t',header=FALSE)

control1 = read.table('file1.tabular',sep='\t',header=FALSE)
control2 = read.table('file2.tabular',sep='\t',header=FALSE)
control3 = read.table('file3.tabular',sep='\t',header=FALSE)

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
jpeg(file="MDSplot_SVA.jpeg", width=5000, height=5000, units="px", res=300)
plotMDS(dge)
dev.off()

#estimate dispersion
disp <- estimateGLMCommonDisp(dge, design)
disp <- estimateGLMTrendedDisp(disp, design)
disp <- estimateGLMTagwiseDisp(disp, design)
plotBCV(disp)

#write plot
jpeg(file="dispersion_plot.jpeg", width=5000, height=5000, units="px", res=300)
plotBCV(disp)
dev.off()


##############################################################

#SVA test#
BiocManager::install("sva")
library(sva)

design0 <- as.data.frame(design)
design0 = model.matrix(~1, data=design0)

#first step
#estimating number of latent factors
n.sv = num.sv(norm_counts,design,method="leek")

#second step
#sva function
svobj = sva(norm_counts,design,design0,n.sv=n.sv)
svobj.df <- data.frame(svobj$sv)

## Source for cleaningP function (see github):
## https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0808-5
cleaningP = function(y, design, svaobj,  P=ncol(design)) {
  X=cbind(design,svaobj$sv)
  Hat=solve(t(X)%*%X)%*%t(X)
  beta=(Hat%*%t(y))
  cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

## Compare PCA with and without SV correction
cleanp = cleaningP(norm_counts,design,svobj)
pca <- prcomp(t(cleanp))
plot(pca)
#autoplot(pca,data=pd,colour='Grupo')

pca0 <- prcomp(t(norm_counts))
plot(pca0)
#autoplot(pca0,data=pd,colour='Grupo')

sv.p3Ascores <- pca$x


######################################################################

##### PCA PLOT #####

#also taken from PCA analysis

jpeg(file="PCA_SVA.jpeg", width=3200, height=3200, units="px", res=300)

## plot the PC
plot(sv.p3Ascores[,1], sv.p3Ascores[,2], xlab="PCA 1", ylab="PCA 2",
     type="p", cex.lab=0.75, cex.axis=0.75, 
     #xlim=c(-200,250), ylim=c(-200,170),
     col=c('deepskyblue','deepskyblue','deepskyblue','darkred','darkred',
           'darkred'),
     main="PCA scores", cex.main=1.2, font.main=1,pch=15)

## apply to plot
text(sv.p3Ascores, colnames(sampleNames), cex=0.5, pos=4, col="black")

legend("topleft", legend=c("Control","COVID"),
       bty="n", xjust = 1, yjust = 1,
       cex=.75, y.intersp=1, col=c('deepskyblue', 'darkred'), pch=20)

##
dev.off()

#######################################################################
##### DIFFERENTIAL EXPRESSION ANALYSIS #####

## Fit the linear model with the surrogate variables included
modSv = cbind(design,svobj.df)

fit <- glmFit(disp, modSv)
lrt <- glmLRT(fit)
topTags(lrt)

#at last ??

# tell edgeR what comparison you want to perform
covidVScontrol = makeContrasts(covid-control, levels=modSv)

# obtain gene symbol
BiocManager::install("org.Hs.eg.db")

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
write.csv(table.covidVScontrol, file="DGE_SVA.csv")
