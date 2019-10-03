########################################################################################
# Functionalized Differential Expression Scatter                                       #
########################################################################################
# Helper functions:                                                                    #
#--------------------------------------------------------------------------------------#
# filter()       <--- optimizeable counts filter                                       #
# nonDegHeat     <--- uses txt output from DE.scatter(nonDEG = TRUE),                  #
#                     maps genes to counts, makes heatmap                              #
# heatClustTxt() <--- Output UHC labels to separate txt files                          #
########################################################################################

library(ggplot2)
library(edgeR)
library(RColorBrewer)
library(ggrepel)
library(DESeq2)
library(gplots) # heatmap.2
library(heatmap3)
library(org.Hs.eg.db)

# Scatter plot of counts, highlights up and down regulated genes (logFC > x or < -x), 
# top 10 up and down regulated genes (normalized counts > 4 and highest logFC score)
#@param
# filteredCounts: read counts file that is filtered as desired
# groupNames: vector of cell types (i.e ['Fibro', 'Fibro', 'iPSC'])
# sizeX:      number of columns for 1st group on X axis
# sizeY:      number of columns for 2nd group on Y axis
# DEGfile:    filepath string of DGEList() to be read, includes logFC in column 1
# keyGenes:   a vector of strings of gene names to be labeled besides top 20 DEGs
# nonDEG:     TRUE = output nonDEGs and their counts to txt file
#@out
# Scatterplot of 
# 'Newtop100.txt': Saves gene names of top 100 upregulated genes to a .txt file for GO analysis
# 'Newbot100.txt': Saves gene names of bottom 100 upregulated genes to a .txt file for GO analysis
# '
DE.scatter <- function(filteredCounts, groupNames, sizeX, sizeY, DEGfile, 
                       logFC = 2, cut = NULL, adjPVal = 0.05, keyGenes =NULL, nonDEG = FALSE) {
  start <- Sys.time()
  deg <- read.delim(DEGfile, row.name = 1)
  X <- filteredCounts[,1:sizeX]
  Y <- filteredCounts[,(sizeX+1):(sizeX+sizeY)]
  # Check if column n=1
  if(sizeX != 1 && sizeY != 1){
    norm_X <- log2(data.frame(Means = rowMeans(X)) + 1)
    norm_Y <- log2(data.frame(Means = rowMeans(Y)) + 1)
  }
  else if(sizeX == 1 && sizeY != 1){
    norm_X <- log2(data.frame(Means = X) + 1)
    norm_Y <- log2(data.frame(Means = rowMeans(Y)) + 1)
  }
  else{
    norm_X <- log2(data.frame(Means = rowMeans(X)) + 1)
    norm_Y <- log2(data.frame(Means = Y) + 1)
  }
  View(cbind(norm_X,norm_Y))
  # Correlation coefficient (using log2(RPM+1))
  cor <- cor.test(norm_X$Means, norm_Y$Means, method = "pearson", conf.level = 0.95)
  cor.coeff <- toString(round(cor$estimate, digits = 4))
  
  # up and downreg gene store
  upreg <- list()
  doreg <- list()
  
  # Color groups for FC > 2, mean counts > 4
  norm_X$color = 'gray'
  norm_X$logFC = ''
  norm_X$adjP = ''
  upCount = 0                               # number of upregulated genes
  doCount = 0                               # number of downregulated genes
  for(i in 1:length(norm_X$Means)){
    FC <- deg[rownames(norm_X[i,]), 1]
    p <- deg[rownames(norm_X[i,]), 5]
    norm_X$logFC[i] <- FC
    norm_X$adjP[i] <- p
    # if(gene logFC < -2  && mean gene counts in B > 4)   
    if(FC < (0 - logFC) && p < adjPVal){
      if(!is.null(cut)){
        if(norm_X[i,1] > cut){
          norm_X[i,2] = 'blue'
          doCount = doCount + 1
          doreg <- rbind(doreg, norm_X[i,])
        }
      }
      else{
        norm_X[i,2] = 'blue'
        doCount = doCount + 1
        doreg <- rbind(doreg, norm_X[i,])
      }
    }
    # else if(gene logFC > 2  && mean gene counts in A > 4)   &&  norm_A[i,1] > 4
    else if(FC > logFC && p < adjPVal){
      if(!is.null(cut)){
        if(norm_Y[i,1] > cut){
          norm_X[i,2] = 'orange'
          upCount = upCount + 1
          upreg <- rbind(upreg, norm_X[i,])
        }
      }
      else{
        norm_X[i,2] = 'orange'
        upCount = upCount + 1
        upreg <- rbind(upreg, norm_X[i,])
      }
    }
  }
  # If nonDEG is TRUE, outpt nonDEG counts to nonDEG.txt
  if(nonDEG){
    print('DEG IF LOOP START')
    grays <- norm_X$color == 'gray'
    nonDEGnames <- rownames(norm_X[grays,])
    print(head(nonDEGnames))
    write(nonDEGnames, file = 'nonDEGnames.txt')
    print('nonDEGs written to: nonDEGnames.txt', quote = FALSE)
  }
  # JOIN iPSC and Fibro means
  norm_X$yMeans <- norm_Y$Means
  norm_X$names <- rownames(norm_X)
  #View(norm_X)
  if(is.null(dim(upreg)) && is.null(dim(doreg))){
    print('No differentially expressed genes', quote = FALSE)
    print(ggplot(norm_X, aes(norm_X$Means, norm_X$yMeans)) + 
            labs(title = paste('correlation coefficient =', cor.coeff), 
                 caption = paste(upCount, 'upreg genes,', doCount, 'downreg genes')) +           # Caption text
            geom_point(color = norm_X$color, size = .5) +                                        # Point formatting
            xlab(paste(groupNames[1], '(log2(RPM+1)')) + ylab(paste(groupNames[2], '(log2(RPM+1)')) + 
            xlim(0,14) + ylim(0,14) +
            theme(aspect.ratio=1, plot.caption = element_text(color = "red", hjust = 0)) +       # Aspect ratio | caption format
            geom_abline(intercept = 0, slope = 1, col = 'black'))
    if(!is.null(keyGenes)){
      print(geom_text_repel(data=subset(norm_X, is.element(names, keyGenes),                     # Label key genes
                                        aes(Means, yMeans, label=names), fontface = 2, point.padding = 1)))
    }
    
  }
  else{
    # Sort up/down reg genes from lowest -> highest logFC
    #View(doreg)
    #View(upreg)
    order.upreg <- upreg[order(as.numeric(upreg$logFC)),]
    order.doreg <- doreg[order(as.numeric(doreg$logFC)),]
    # Write gene names to .txt file for GO analysis later
    write(rownames(tail(order.upreg, 100)), file = 'NewUp100.txt')
    write(rownames(head(order.doreg, 100)), file = 'NewDown100.txt')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', quote = FALSE)
    print('Top 100 upregulated and downregulated gene names out put to:', quote = FALSE)
    print('NewUp100.txt     NewDown100.txt', quote = FALSE)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', quote = FALSE)
    top_reg <- c(rownames(tail(order.upreg, 10)), rownames(head(order.doreg, 10)))
    if(!is.null(keyGenes)){
      top_reg <- append(top_reg, keyGenes)
    }
    # Plot configs
    print(ggplot(norm_X, aes(norm_X$Means, norm_X$yMeans)) + 
            labs(title = paste('correlation coefficient =', cor.coeff), 
                 caption = paste(upCount, 'upreg genes,', doCount, 'downreg genes')) +           # Caption text
            geom_point(color = norm_X$color, size = .5) +                                        # Point formatting
            xlab(paste(groupNames[1], '(log2(RPM+1)')) + ylab(paste(groupNames[2], '(log2(RPM+1)')) + 
            xlim(0,14) + ylim(0,14) +                                                            #TODO
            theme(aspect.ratio=1, plot.caption = element_text(color = "red", hjust = 0)) +       # Aspect ratio | caption format
            geom_abline(intercept = 0, slope = 1, col = 'black') +
            geom_text_repel(data=subset(norm_X, is.element(names, top_reg)),           # Label top genes
                            aes(Means, yMeans, label=names), fontface = 1, point.padding = 1))
  }
  end <- Sys.time()
  print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', quote = FALSE)
  print(paste('p-value:', adjPVal), quote = FALSE)
  print(paste('Count cutoff:', cut), quote = FALSE)
  print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', quote = FALSE)
  print(end-start)
}

#  Counts data cleanup helper function
filter <- function(seq.data, cpmThresh = 0.5, minCount = NULL){
  seq.data <- data.frame(seq.data)
  if(is.null(minCount)){
    min <- dim(seq.data)[2]/2
  }
  else{min <- minCount}
  seq_cpm <- cpm(seq.data)
  thresh <- seq_cpm > cpmThresh
  keep <- rowSums(thresh) >= (min)
  keep <- seq_cpm[keep,]
  #print(head(rownames(keep)))
  #print(head(colnames(keep)))
  
  return(keep)
}

# Helper function for if DE.scatter() was run with nonDEG = TRUE
# collects nonDEG counts and makes heatmap
nonDegHeat <- function(nonDEG.txt, counts, title){
  nonDEGs <- read.table(nonDEG.txt)
  nonD <- unlist(nonDEGs, use.names = FALSE)
  issaNondeg <- rownames(counts) %in% nonD
  nonDegCounts <- counts[issaNondeg,]
  # Make heatmap
  y <- DGEList(nonDegCounts)
  y <- calcNormFactors(y)
  logcounts <- cpm(y,log=TRUE)
  
  quantile.range <- quantile(logcounts, probs = seq(0, 1, 0.01))
  palette.breaks <- seq(quantile.range["10%"], quantile.range["90%"], 0.1)
  pal  <- colorRampPalette(c("red",'white', "blue"))(length(palette.breaks) - 1)
  more <- colorRampPalette(pal)
  rowLabs <- c(rep('', length(nonD)))
  title <- paste('Non-DEGs Heatmap:', length(nonD), 'genes')
  heatmap3(logcounts)
  title <- paste(title)
  heatmap3(logcounts, main = title)
}

# PCA --> Loadings Scatter --> Heatmap --> clustered gene UHC
# Helper function to allocate cluster genes in separate txt files
heatClustTxt <- function(folderPath, heatmap, numClust){
  heatUHC <- as.hclust(heatmap$rowDendrogram)
  clusters <- data.frame(cutree(heatUHC, numClust))
  clustNames <- cbind(clusters, rownames(clusters))
  View(clustNames)
  for (i in 1:numClust){
    filePath <- paste(folderPath, 'clust', i, '.txt', sep = '')
    condition <- clusters == i
    clust <- data.frame(condition[condition,])
    write.table(rownames(clust), file = filePath, row.names = FALSE, col.names = FALSE, quote = FALSE)
    print(paste('Cluster', i ,'gene names saved to:', filePath))
  }
}

# PCA --> Loadings Scatter --> Heatmap --> clustered gene UHC --> ENTREZ IDs for topGene
# Helper function for converting each cluster txt from heatClustTxt() to entrez txt
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
