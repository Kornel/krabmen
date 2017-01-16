library(readr)
library(dplyr)
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

