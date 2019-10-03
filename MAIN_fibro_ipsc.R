# MAIN: Rhino fibro and ipsc
## Paper has average across sample counts, log2(RPM+1)
## The +1 is for the genes that have logcount = 0
## >4-fold difference upper and lower bound
library(pheatmap)
library(RColorBrewer)
library(ggdendro)
library(ggplot2)
library(edgeR)
library(ggrepel)
library(rgl)
library(gplots) # heatmap.2
library(factoextra)
library(dendextend)
library(viridis)
library(animation)
library(magick)

seqdata <- read.delim("~/ConGen/counts/all_Fibro_and_iPSCs.txt", row.names=1)
fibro.ipsc.counts <- filter(seqdata)
head(fibro.ipsc.counts)
# Generate DEG list with logFC
degOut <- '~/ConGen/logFC/fibro(-)_ipsc(+).txt'

fibro <- fibro.ipsc.counts[,1:7]
ipsc <- fibro.ipsc.counts[,8:14]

y <- DGEList(fibro.ipsc.counts)
y <- calcNormFactors(y) 

group <- c(rep("fibro",7), rep("ipsc", 7))
group <- as.factor (group)
design <- model.matrix(~ 0 + group)

v <- voom(y, design, plot = TRUE)
fit <- lmFit(v)

cont.matrix <- makeContrasts(groupipsc - groupfibro,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
cont.matrix
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)
limma.res <- topTable(fit.cont,sort.by="p",n="Inf")
limma.res <- topTable(fit.cont,n="Inf")
#write.table(limma.res,file=fileOut,sep="\t")

# Pairwise DEG Scatterplot
DE.scatter(filter(fibro.ipsc.counts), c("Fibroblasts", "iPSCs"), 7, 7, 
           DEGfile = '~/ConGen/logFC/fibro(-)_ipsc(+).txt')
# keyGenes = c('TFAP2A', 'FBN1', 'MTG1'), nonDEG = TRUE)

# Heatmap of nonDEGs from DEG scatterplot
nonDegHeat("~/ConGen/nonDEGnames.txt", fibro.ipsc.counts)


# UHC Samples
hc <- hclust(dist(t(fibro.ipsc.counts)), 'average')
color.uhc <- hc %>%
  color_branches(k = 2) %>%
  set("branches_lwd", c(2,1,2)) %>%
  set("branches_lty", c(1,2,1))
par(mar=c(3,1,3,8))
plot(color.uhc, horiz = TRUE)
par(mar=c(4,4,4,4))


# 2D PCA
y <- DGEList(fibro.ipsc.counts)
y <- calcNormFactors(y)                    # edgeR normalization
cell.type <- c(rep("Fibroblasts",7), rep("iPSCs",7))
norm <- cpm(fibro.ipsc.counts)
pca <- prcomp(t(norm))                      # prcomp requires transformed counts dataframe
fviz_pca_ind(pca, col.ind = cell.type, repel = TRUE, invisible = 'quali')

# 3D PCA
type = c(rep('blue', 7), rep('orange', 7))
plot3d(pca$x, col = type, type = 's', size = 1,
       xlab = 'PC1', ylab = 'PC2', zlab = 'PC3') +
       grid3d('z', at = NULL, col = "gray", lwd = 1, lty = 1, n = 5)                 # add grid
legend3d("topright", legend = c('Fibroblast', 'iPSC'), 
         pch = 16, col = c('blue', 'orange'), cex=2, inset=c(0.02))
text3d(pca$x, texts = colnames(fibro.ipsc.counts))

notrun <- "
tempdir()
movie3d(spin3d(axis = c(0, 1, 0), rpm = 3), convert = TRUE, duration = 20)
# Capture snapshot to filename (Will overwrite filename if left unchanged!)
snapshot3d(filename = 'pictures/3dPCA(1)_labels.png', fmt = 'png')
"

# PCA Loadings Scatterplot (using scale and center params? Would look more uniform circle)
# Notes:
# PCAdist is a helper function in the FUNC_scatter.R file
rotation <- data.frame(pca$rotation[,1:2])
distance <- function(x, y) sqrt(x^2 + y^2)
rotation <- rotation * 1000
head(rotation)
rotation$color <- 'gray'
radius <- 2
pLims <- 8     # Scatterplot window limits
count <- 0     # Tally number of significant genes
for (i in 1:nrow(rotation)){
  if (distance(rotation[i,1], rotation[i,2]) > radius) {
    rotation$color[i] <- 'red'
    count <- count+1
  }
}

#TODO get top 5 genes from 4 quadrants
rotation <- rotation[order(as.numeric(rotation$PC1)),]

plot(rotation[,1:2], pch=20, cex = .5,  
     col = color, xlim=c(-pLims, pLims), ylim=c(-pLims, pLims)) +
     mtext(paste('rPC1:2 > SD2 |', count, 'genes'), side = 3) + 
     grid(col = "lightgray", lty = 'dashed', nx =20, ny = 20)
text(rotation[,1:2], labels = rotation$names, pos = 1)

# Heatmap from loading scatter genes
# pretty colorbrewer palettes: viridis, magma, plasma, inferno
# At a glance, magma looks like most universally colorblind friendly
# For clear instructions on how to adjust gradient: 
# https://ryjohnson09.netlify.com/post/how-to-make-a-heatmap-in-r/
selected <- rotation$color == 'red'
gene.select <- rownames(rotation[selected,])
gene.select.data <- fibro.ipsc.counts[gene.select,]
head(gene.select.data)
y <- DGEList(gene.select.data)
y <- calcNormFactors(y)
logcounts <- cpm(y,log=TRUE)
rpm1 <- logcounts +1
numRows <- dim(gene.select.data)[1]
title <- 'NWR Fibroblasts and iPSCs: 
          PCA loading score significant genes'
sampleCol <- c(rep('steel blue', 7), rep('orange',7))
rowLabs <- c(rep('', dim(rpm1)[1]))

quantile.range <- quantile(logcounts, probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["10%"], quantile.range["90%"], 0.1)
pal <- colorRampPalette(c('blue','white', 'red'))

heat <- heatmap.2(rpm1,col=magma, breaks = c(1, 5:11, seq(12,max(rpm1), 5), 20),
                  main=title, srtCol = 45, trace="none",
                  keysize = 1, key.title = paste('log2(RPM+1)'), 
                  labRow = rowLabs, cexRow = 1,
                  ColSideColors = sampleCol, colsep = 7, sepcolor = 'white', sepwidth = 0.05)

# Heatmap genes UHC
heat.clust.genes <- rownames(logcounts)[heat$rowInd]
head(heat.clust.genes)
#file <- 'heat.uhc.genes.fibro.ipsc.txt'
write.table(heat.clust.genes, file = file, row.names = FALSE, col.names = FALSE, quote = FALSE)
print(paste('Heatmap UHC gene names output to:', file))

heatUHC <- as.hclust(heat$rowDendrogram)
View(heatUHC)

# output clusters of gene names in separate txt files
heatClustTxt('~/ConGen/DAVID/DAVID_fibro_ipsc/', heat, numClust = 5)
# convert to Entrez IDs for topGene
library(org.Hs.eg.db)
clustToEntrez <- function(folderPath, numClust){
  hs <- org.Hs.eg.db
  for(i in 1:numClust){
    dest <- paste(folderPath, '/entrez', i,'.txt',sep='')
    cluster <- read.table(paste(folderPath, '/clust', i, '.txt', sep = ''), quote="\"", comment.char="")
    vector <- as.vector(cluster[,'V1'])
    entrez <- select(hs, keys = vector, columns = "ENTREZID", keytype = "SYMBOL")
    write.table(entrez[,2], file = dest, row.names = FALSE, col.names = FALSE, quote = FALSE)
    print(paste('Cluster', i, 'Entrez IDs output to:', dest))
  }
}
clustToEntrez('~/ConGen/DAVID/DAVID_fibro_ipsc', 5)

# Heatmap dendrogram of gene clusters
color.uhc <- heat.uhc %>%
  color_branches(k = 5) %>%
  set("branches_lwd", c(2,1,2)) %>%
  set("branches_lty", c(1,2,1))
dend <- color_labels(color.uhc, k = 5, par = 2)
par(mar=c(0,0,0,6))
plot(dend, horiz = TRUE)
par(mar=c(4,4,4,4))

# GO analysis using topGO
library(topGO)
library(ALL)
data(ALL)
data(geneList)
library(topGO)
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)

DEGfile = read.delim('~/ConGen/logFC/fibro(-)_ipsc(+).txt', row.names = 1)
DEG <- as.vector(DEGfile[,1])
names(DEG) <- rownames(DEGfile)

genes <- rownames(gene.select.data[1:50,])
ficounts <- genes
typeof(ficounts)
ficounts <- as.vector(ficounts)
sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = DEG, geneSel = topDiffGenes, 
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)
