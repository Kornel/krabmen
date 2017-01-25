library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
file <- '../data/methylation/BRCA.methylation.27k.450k.txt.gz'

rawdata <- read_delim(file, delim = '\t')

#Only complete cases
data <- rawdata[complete.cases(rawdata),]

#Short colnames
colnames(data) <- substr(colnames(data), 0, 12)

#Subtypes data for clusters
subtypes <- read_delim('../KRAB_and_TCGA_Subtypes/Data v3 2016.06.06/BRCA.csv', delim = ';')
subtypes <- subtypes[, c('sampleID', 'methylation.Clusters')]
subtypes <- subtypes[complete.cases(subtypes),]

#Transpose probes to columns
probes <- data$id
samples <- colnames(data)

#Remove probe names
data <- data[,-1]

#Transpose
tdata <- as.data.frame(t(data))

#Join with subtypes to get methylation.Cluster
tdata$sampleID <- rownames(tdata)
jdata <- tdata %>% inner_join(subtypes, by = c('sampleID' = 'sampleID'))

colnames(jdata) <- c(probes, c('sampleID', 'cluster'))

jdata$cluster <- as.factor(jdata$cluster)

# Density
library(ggplot2)
ggplot(jdata) + geom_boxplot(aes(x=cluster, y=cg00000292), alpha = 0.4)
ggplot(jdata) + geom_boxplot(aes(x=cluster, y=cg00002426), alpha = 0.4)

# Sample data
features <- 1000
observations <- 935
head <- jdata[complete.cases(jdata),]
rdu <- function(n, min, max) sample(min:max, n, replace = F)
head <- head[rdu(observations, 1, dim(jdata)[1]),] 
head$sampleID <- NULL
idx <- length(colnames(head))
head <- head[, c(rdu(n = features, min = 1, max = idx - 1), idx)]
idx <- length(colnames(head))

# Andrews plot on sample
andrews(head, clr = idx)

# PCA on sample
log.head <- log(head[, 1:(idx - 1)])
cluster <- head[, idx]

pca <- prcomp(log.head, center = T, scale. = T)
plot(pca, type = 'l')

pca.df <- as.data.frame(as.matrix(log.head) %*% pca$rotation)
pca.df$cluster <- cluster
ggplot(pca.df) + geom_point(aes(x = PC1, y = PC2, colour = cluster))

# PCA on all
log.jdata <- log(jdata[, 1:(dim(jdata)[2] - 2)])
jdata.cluster <- jdata[, dim(jdata)[2]]

pca <- prcomp(log.jdata, center = T, scale. = T)
plot(pca, type = 'l')

pca.df <- as.data.frame(as.matrix(log.jdata) %*% pca$rotation)
pca.df$cluster <- jdata.cluster
ggplot(pca.df) + geom_point(aes(x = PC1, y = PC2, colour = cluster))

(sums <- cumsum(pca$sdev^2 / sum(pca$sdev^2)))

library(car)
scatter3d(pca.df$PC3, pca.df$PC2, pca.df$PC1, groups = pca.df$cluster, surface = F)


# Clusters
library(cluster)
# KMeans
kmAll <- kmeans(pca.df, centers = 3)
km2d <- kmeans(pca.df[,c(1,2)], centers = 3)

# KMedoids
pamAll <- pam(pca.df[, c(1, 2)], k = 3)

# GMM
library(mclust)
mclust <- Mclust(pca.df)
plot(xyMclust)

#3d scatterplot with clusters
scatter3d(pca.df$PC3, pca.df$PC2, pca.df$PC1, groups = as.factor(km$cluster), surface = F)

#2d scatterplot with clusters
ggplot(pca.df) + geom_point(aes(x = PC1, y = PC2, colour = as.factor(pamAll$clustering)))

# Save clusters assigned to probes
pca.clusters <- data.frame(
  sampleID = jdata$sampleID,
  clustser2PCs = km2d$cluster
)
