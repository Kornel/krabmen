library(readr)

source('./data-utils.R')
source('./utils.R')

files <- select.tumor.files(patient.freqs.path = '../../results/stats/normalized/stats-rsem-normalized-1-vs-11.csv',
                            files.path = '../../download/mrna-rsem-normalized-2015-11-01/')

file <- files[1]

for (file in files) {
  tumor.name <- filename.to.tumor(file)
  
  print(tumor.name)
  rawdata <- read_delim(file, delim = '\t')  
  rawdata$HybRefShort <- hyb.ref.short(rawdata)
  rawdata$`Hybridization REF` <- NULL
  gene.table <- get.gene.table(rawdata, tumor.name)$gene.table
  
  gene.table <- t(gene.table)
  r <- rownames(gene.table)
  rownames(gene.table) <- NULL
  gene.table <- as.data.frame(gene.table)
  gene.table$type <- r
  
  gene.table <- gene.table %>% filter(type != 'other')

  save(gene.table, file = paste0('../../../boxplots-krabmen/resources/', tumor.name, '.Rdata'))  
}

