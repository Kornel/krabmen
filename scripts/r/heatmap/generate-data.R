missing.genes <- data.frame(tumor = character(0), gene.id = character(0))

heatmap.data <- data.frame(log2median = numeric(0), 
                           log2mean = numeric(0), 
                           Gene.Id = character(0), 
                           tumor.name = character(0))



rawdata <- read.csv('../../results/bionmialTest/full-table-pvalues-long.csv')

selected <- subset(rawdata, select = c('id', 'tumor', 'medianHealthy', 'medianTumor', 'log2FoldChange'))
colnames(selected) <- c('Gene.Id', 'tumor.name', 'medianHealthy', 'medianTumor', 'log2FoldChangeMean')

selected$log2FoldChangeMedian <- log2(selected$medianTumor / selected$medianHealthy)

selected <- subset(selected, select = c('Gene.Id', 'tumor.name', 'log2FoldChangeMedian', 'log2FoldChangeMean'))

write.csv(selected, file = sprintf('../../results/heatmap/heatmap-data.csv'), row.names = F)
