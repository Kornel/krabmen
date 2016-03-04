library(DESeq)
library(dplyr)
library(readr)
library(reshape2)
library(data.table)
library(pheatmap)


genes <- read.csv('../../data/KRAB ZNF gene master list.csv', sep = '\t')
genes$Gene.ID <- as.character(genes$Gene.ID)

source('utils.R')

patient.freqs.path <- '../../results/stats/raw/stats-rsem-1-vs-11.csv'
patient.freqs <- read.csv(patient.freqs.path) %>% filter(normal >= 9)

files.path <- '../../download/mrna-rsem-2015-11-01'

files <- list.files(path = files.path, pattern = '*-clean.txt', recursive = T)

filtered.files <- c()

for (tumor in patient.freqs$tumor.name) {
  contains <- lapply(files, function(x) grepl(tumor, x))
  file <- files[contains == T]
  filtered.files <- c(file, filtered.files)
}

genes <- read.csv('../../data/KRAB ZNF gene master list.csv', sep = '\t')
genes$Gene.ID <- as.character(genes$Gene.ID)

source('utils.R')


for (file in filtered.files) {
  tumor.name <- sub('gdac.broadinstitute.org_(\\w*)\\..*', '\\1', file)
  
  print(tumor.name)
  
  x <- read_delim(sprintf('%s/%s', files.path, file), delim = '\t')
  
  x2 <- x[,seq(2, length(x), 3)]
  
  x2$HybRefShort <- sub('\\|.*', '', x$`Hybridization REF`)
  
  x2 <- genes %>% inner_join(x2, by = c('Gene.ID' = 'HybRefShort'))
  x2$HybRefShort <- NULL
  x2$Ensembl.ID <- NULL
  x2$`Hybridization REF` <- NULL
  
  rownames(x2) <- x2$Gene.ID
  
  x2$Gene.ID <- NULL
  
  x3 <- as.data.frame(sapply(x2, as.integer))
  rownames(x3) <- rownames(x2)
  
  condition <- factor(barcode.to.type(colnames(x3)))
  
  count2 <- x3[,condition != 'other']
  cond2 <- factor(barcode.to.type(colnames(count2)))
  cds <- newCountDataSet(count2, cond2)
  cds <- estimateSizeFactors(cds)
  
  cds <- estimateDispersions(cds)
  #plotDispEsts(cds)
  
  res <- nbinomTest(cds, 'healthy', 'tumor')
  write.csv(res, file = sprintf('../../results/bionmialTest/%s.csv', tumor.name), row.names = F)
}

files <- list.files(path = '../../results/bionmialTest', pattern = '*.csv', recursive = T, full.names = T)

all <- data.frame()



for (file in files) {
  tumor.name <- sub('.*/(\\w*)\\.csv', '\\1', file)
  print(tumor.name)
  
  raw <- read.csv(file)
  raw$tumor <- tumor.name
  all <- rbind(all, raw)
}

values <- c('baseMean','baseMeanA', 'baseMeanB', 'foldChange', 'log2FoldChange', 'pval', 'padj')
all.2 <- as.data.frame(dcast(setDT(all), id ~ tumor, value.var = values))

cnames <- colnames(all.2)[-1]
order <- c(1, order(gsub(".*_(.*)", "\\1", cnames)) + 1)

cnames <- c('id', cnames[order(gsub(".*_(.*)", "\\1", cnames))])

sorted <- all.2[,cnames]

write.csv(all, file = '../../results/bionmialTest/heatmap-data.csv', row.names = F)
write.csv(sorted, file = '../../results/bionmialTest/full-table-pvalues.csv', row.names = F)

heatmap.data <- read.csv('../../results/bionmialTest/heatmap-data.csv')

h <- heatmap.data[, c('tumor', 'id', 'log2FoldChange')]  

h.mean <- dcast(h, tumor ~ id, value.var = 'log2FoldChange')
rownames(h.mean) <- h.mean$tumor.name
h.mean$tumor.name <- NULL
h.mean <- t(h.mean)
pheatmap(h.mean, 
         show_colnames = T, 
         cluster_cols = T, 
         cluster_rows = T, 
         fontsize_row = 2,
         #   filename = sprintf('%s/heatmap-mean.png', results.dir),
         width = 10,
         height = 20)

