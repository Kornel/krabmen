library(readr)

source('./data-utils.R')
source('./utils.R')
source('./rsem-normalized/stats.R')

tumor.files <- select.tumor.files(
  patient.freqs.path = '../../results/stats/normalized/stats-rsem-normalized-1-vs-11.csv',
  files.path = '../../download/mrna-rsem-normalized-2015-11-01/'
)

file <- tumor.files[1]

for (file in tumor.files) {
  tumor.name <- filename.to.tumor(file)
  
  print(tumor.name)
  
  rawdata <- read_delim(file, delim = '\t')
  #Remove second header
  data <- rawdata[-1, ]
  #All except first column - Hybridization REF
  data[, -1] <- as.data.frame(sapply(data[, -1], as.numeric))
  
  data$HybRefShort <- hyb.ref.short(data)
  data$`Hybridization REF` <- NULL
  filtered.counts <- get.gene.table(data, tumor.name)$gene.table
  
  counts <- filtered.counts
  counts <- counts[, colnames(counts) != 'other']

  res <- get.full.stats(counts)
  
  write.csv(res, file = sprintf('../../results/rsem-normalized/partial/%s.csv', tumor.name), row.names = F)
}
