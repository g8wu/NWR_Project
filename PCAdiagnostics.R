# PCA diagnostics
# All of samples have very high variance except 1
# which means 14 samples can explain more htan 90% of variance
# You can reduce dimensionality from 15 to 14 and lose 10% variance
summary(scaled.pca)
dim(scaled.pca$rotation)

normalize <- function(vector){
  out <-c()
  min <- min(vector)
  max <- max(vector)
  if(min < 0){ vector <- abs(min) + vector}
  out <- (vector-min)/(max-min)
  return(t(out))
}
test.pca <- scaled.pca$rotation[1,]
normalize(test.pca)

norm.pca <- data.frame()
for(i in 1:dim(scaled.pca$rotation)[1]){
  norm.vect <- normalize(scaled.pca$rotation[i,])
  norm.pca <- rbind(norm.pca, norm.vect)
}
dim(norm.pca)
rownames(norm.pca)<- rownames(scaled.pca$rotation)
View(norm.pca)

norm.cat.pca <- norm.pca[1:2]
norm.cat.pca$color <- 'gray'
for (i in 1:length(norm.cat.pca[,1])){
  if (dist(norm.cat.pca[i,1], norm.cat.pca[i,2]) > radius) {
    norm.cat.pca[i, 3] <- 'red'
    count <- count+1
  }
}
caption <- paste('rPC1_2 > SD', radius, ' | ', count, ' genes', sep = '')
plot(norm.cat.pca[,1:2], pch=20, cex = .5, main="Normalized Loading Scores", 
     col = rotation[,3])
mtext(caption, side = 3)
draw.circle(0,0,radius,border="blue", col = NA)





screeplot(scaled.pca, type = "l", npcs = 15, main = "Screeplot of PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
# You can explain almost 80% of variance with first 6 components
cumpro <- cumsum(scaled.pca$sdev^2 / sum(scaled.pca$sdev^2))
plot(cumpro, xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.88759, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)