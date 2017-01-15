library(readr)
library(dplyr)
library(reshape2)
file <- 'data/out.txt'
#file <- 'data/head.txt'

rawdata <- read_delim(file, delim = '\t')
#Remove second header
data <- rawdata[-1, ]

#Beta values are in every 4th column, 
#first column is the name of the probe
cols <- c(1, seq(2, dim(data)[2], 4))
data <- data[, cols]

#Convert Beta values to numeric
data[, -1] <- as.data.frame(sapply(data[, -1], as.numeric))

#Only complete cases
data <- data[complete.cases(data),]

#Short colnames
colnames(data) <- substr(colnames(data), 0, 12)

#Subtypes data for clusters
subtypes <- read_delim('../KRAB_and_TCGA_Subtypes/Data v3 2016.06.06/BRCA.csv', delim = ';')
subtypes <- subtypes[, c('sampleID', 'methylation.Clusters')]
subtypes <- subtypes[complete.cases(subtypes),]

#Transpose probes to columns
probes <- data$Hybridizatio
samples <- colnames(data)

#Remove probe names
data <- data[,-1]

#Transpose
tdata <- as.data.frame(t(data))

#Join with subtypes to get methylation.Cluster
tdata$sampleID <- rownames(tdata)
jdata <- tdata %>% inner_join(subtypes, by = c('sampleID' = 'sampleID'))

colnames(jdata) <- c(probes, c('sampleID', 'cluster'))

#Anova

.append.stat <- function(valueCol, discriminator, stat, name, df) {
  stats <- tapply(valueCol, discriminator, stat)
  names(stats) <- paste0(name, names(stats))
  cbind(df, t(as.data.frame(stats)))
}

result <- data.frame()

for (probe in probes) {
  probeCol <- jdata[, probe]
  clusterCol <- as.character(jdata$cluster)
  # Anova
  a <- aov(probeCol ~ clusterCol)
  
  #TukeyHSD
  t.hsd <- as.data.frame(TukeyHSD(a)$`clusterCol`)
  
  # Convert TukeyHSD result to dataframe
  t.hsd$type <- rownames(t.hsd)
  t.hsd$probe <- probe
  t.hds.df <- dcast(t.hsd, probe ~ type, value.var = c('p adj'))
  colnames(t.hds.df)[-1] <- paste0('p-adj-', colnames(t.hds.df)[-1])
  
  # Append stats: mean, SEM, Q1, Q3
  t.hds.df <- .append.stat(probeCol, clusterCol, mean, 'mean-', t.hds.df)
  t.hds.df <- .append.stat(probeCol, clusterCol, function(x) sqrt(var(x) / length(x)), 'SEM-', t.hds.df)
  t.hds.df <- .append.stat(probeCol, clusterCol, function(x) quantile(x, 1/4), 'Q1-', t.hds.df)
  t.hds.df <- .append.stat(probeCol, clusterCol, function(x) quantile(x, 3/4), 'Q3-', t.hds.df)
  t.hds.df <- .append.stat(probeCol, clusterCol, length, 'count-', t.hds.df)

  result <- rbind(result, t.hds.df)
}
