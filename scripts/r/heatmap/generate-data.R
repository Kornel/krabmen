# Prepare heatmap data from deseq normalized data

rawdata <- read.csv('../../results/bionmialTest/full-table-pvalues-long.csv')

rawdata <- rawdata[rawdata$tumor != 'KIPAN',]

selected <- subset(rawdata, select = c('id', 'tumor', 'medianHealthy', 'medianTumor', 'log2FoldChange'))
colnames(selected) <- c('Gene.Id', 'tumor.name', 'medianHealthy', 'medianTumor', 'log2FoldChangeMean')

selected$log2FoldChangeMedian <- log2(selected$medianTumor / selected$medianHealthy)

selected <- subset(selected, select = c('Gene.Id', 'tumor.name', 'log2FoldChangeMedian', 'log2FoldChangeMean'))

selected <- do.call(data.frame, lapply(selected, function(x) replace(x, is.infinite(x), NA)))
write.csv(selected, file = sprintf('heatmap/data/heatmap-data-deseq-normalized.csv'), row.names = F)

# Prepare heatmap data from rsem normalized data

rawdata <- read.csv('../../results/rsem-normalized/full-table-long.csv')
rawdata <- rawdata[rawdata$tumor != 'KIPAN',]

selected <- subset(rawdata, select = c('Gene.ID', 'tumor', 'log2FoldChangeMean', 'log2FoldChangeMedian'))
colnames(selected)[1] <- 'Gene.Id'
colnames(selected)[2] <- 'tumor.name'

selected <- do.call(data.frame, lapply(selected, function(x) replace(x, is.infinite(x), NA)))

write.csv(selected, file = sprintf('heatmap/data/heatmap-data-rsem-normalized.csv'), row.names = F)
