library(DESeq)
library(dplyr)
library(readr)

source('./data-utils.R')
source('./utils.R')
source('./binomtest/stats.R')

tumor.files <- select.tumor.files()

file <- tumor.files[1]

for (file in tumor.files) {
  tumor.name <- filename.to.tumor(file)
  
  print(tumor.name)
  
  rawdata <- read_delim(file, delim = '\t')
  
  raw.counts <- raw.counts.only(rawdata)
  
  filtered.counts <- get.gene.table(raw.counts, tumor.name)$gene.table

  counts <- as.data.frame(sapply(filtered.counts, as.integer))
  rownames(counts) <- rownames(filtered.counts)
  
  counts <- counts[, colnames(counts) != 'other']
  conditions <- factor(gsub('(.*)\\..*', '\\1', colnames(counts)))
  
  cds <- newCountDataSet(counts, conditions)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  res <- nbinomTest(cds, 'healthy', 'tumor')
  
  res <- res %>% inner_join(get.stats(counts), by = c('id' = 'Gene.ID'))
  res <- plyr::rename(res, c('baseMeanA' = 'baseMeanHealthy', 'baseMeanB' = 'baseMeanTumor'))
  
  write.csv(res, file = sprintf('../../results/bionmialTest/partial/%s.csv', tumor.name), row.names = F)
}



colnames(counts) factor(gsub('(.*)\\..*', '\\1', colnames(counts)))

counts

tumor <- counts[1, gsub('(.*)\\..*', '\\1', colnames(counts)) == 'tumor'] 
healthy <- counts[1, gsub('(.*)\\..*', '\\1', colnames(counts)) == 'healthy'] 

t <- unlist(unname(as.list(tumor)))
h <- unlist(unname(as.list(healthy)))

type <- c(rep('t', length(t)), rep('h', length(h)))
boxplot(c(t,h) ~ type)


hist(h)
