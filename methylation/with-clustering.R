library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(car)
library(RColorBrewer)
library(mclust)

# Analyze p-values

results.A <- read.csv('results/BRCA-with-probe-mapping.csv', sep = ',')

pvalues.A <- results.A[, c('probe', colnames(results.A)[grepl('p.adj.*', colnames(results.A))])]

pvalues.A <- melt(pvalues.A, id = c('probe'))

unusual.pvals <- c('p.adj.3.2', 'p.adj.5.3', 'p.adj.4.3', 'p.adj.3.1')
pvalues.A.f <- pvalues.A %>% filter(!(variable %in% unusual.pvals))
pvalues.A.f2 <- pvalues.A %>% filter(variable %in% unusual.pvals)

ggplot(pvalues.A.f) + 
  geom_density(aes(value, fill = variable), colour = F, alpha = 0.5)

ggplot(pvalues.A.f2) + 
  geom_histogram(aes(value, fill = variable), position = 'dodge', bins = 30)


# PCA, Cluster

file <- '../data/methylation/BRCA.methylation.27k.450k.txt.gz'

rawdata <- read_delim(file, delim = '\t')

#Only complete cases
data <- rawdata[complete.cases(rawdata),]

#Transpose probes to columns
probes <- data$id

#Remove probe names
data <- data[,-1]

#Transpose
tdata <- as.data.frame(t(data))
colnames(tdata) <- probes

#Log skewed data
log.data <- log(tdata)

pca <- prcomp(log.data, center = T, scale. = T)
plot(pca, type = 'l')

dim(log.data)
dim(pca$rotation)

# pca.df <- as.data.frame(as.matrix(log.data) %*% pca$rotation)
pca.df <- as.data.frame(pca$x)

(sums <- cumsum(pca$sdev^2 / sum(pca$sdev^2)))

# Take a look in 2d
ggplot(pca.df) + geom_point(aes(x = PC1, y = PC2))

# Take a look in 3d
scatter3d(pca.df$PC3, pca.df$PC2, pca.df$PC1, surface = F)

# Clustering
# Determine number of clusters

wss <- (nrow(pca.df)-1) * sum(apply(pca.df, 2, var))
max.clusters <- 200
for (i in 2:max.clusters) {
  if (i %% 10 == 0) print(i / max.clusters)
  wss[i] <- sum(kmeans(pca.df, centers = i)$withinss)
}
plot(1:max.clusters, wss, type = 'b', xlab = 'Number of Clusters', ylab = 'Within groups sum of squares')

km <- kmeans(pca.df, centers = 50)

clusters <- as.factor(km$cluster)

# Take a look in 2d
ggplot(pca.df) + geom_point(aes(x = PC1, y = PC2, colour = clusters))

# Take a look in 3d
colors <- colorRampPalette(brewer.pal(9, 'Set1'))(50)
scatter3d(pca.df$PC3, pca.df$PC2, pca.df$PC1, group = clusters, surface.col = colors, point.col = colors, surface = F)

# Mclust
clust <- Mclust(pca.df, G = 50)
# summary(clust, parameters = TRUE)
ggplot(pca.df) + geom_point(aes(x = PC1, y = PC2, colour = as.factor(clust$classification)))
ggplot(pca.df) + geom_point(aes(x = PC1, y = PC3, colour = as.factor(clust$classification)))
ggplot(pca.df) + geom_point(aes(x = PC2, y = PC3, colour = as.factor(clust$classification)))
scatter3d(pca.df$PC3, pca.df$PC2, pca.df$PC1, group = as.factor(clust$classification), surface.col = colors, point.col = colors, surface = F)


sampleIDs <- rownames(pca.df)
pca.test <- tdata
pca.test$cluster <- as.factor(clust$classification)

result <- data.frame()

components <- head(colnames(pca.test), -1)

n <- 1
for (cmp in components) {
  
  print(n / length(components) * 100)
  n <- n + 1
  compCol <- pca.test[, cmp]
  clusterCol <- as.character(pca.test$cluster)
  # Anova
  a <- aov(compCol ~ clusterCol)
  
  #TukeyHSD
  t.hsd <- as.data.frame(TukeyHSD(a)$`clusterCol`)
  
  # Convert TukeyHSD result to dataframe
  t.hsd$type <- rownames(t.hsd)
  t.hsd$component <- cmp
  t.hds.df <- dcast(t.hsd, component ~ type, value.var = c('p adj'))
  colnames(t.hds.df)[-1] <- paste0('p-adj-', colnames(t.hds.df)[-1])
  
  result <- rbind(result, t.hds.df)
}

write.table(result, 'results/pca-clusters.csv', row.names = F, sep = ',')

results.B <- read.csv('results/pca-clusters.csv', sep = ',')

pvalues.B <- results.B[, c('component', colnames(results.B)[grepl('p.adj.*', colnames(results.B))])]

pvalues.B <- melt(pvalues.B, id = c('component'))

pvalues.B %>% filter(variable %in% unique(pvalues.B$variable)[100:120]) %>% ggplot() + 
  geom_density(aes(value, fill = variable), colour = F, alpha = 0.5)

ggplot(pvalues.A.f2) + 
  geom_histogram(aes(value, fill = variable), position = 'dodge', bins = 30)

results.B$
