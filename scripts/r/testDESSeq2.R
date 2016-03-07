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


file <- filtered.files[1]

for (file in filtered.files) {
  tumor.name <- sub('gdac.broadinstitute.org_(\\w*)\\..*', '\\1', file)
  
  print(tumor.name)
  
  rawdata <- read_delim(sprintf('%s/%s', files.path, file), delim = '\t')
  
  raw.counts <- rawdata[,seq(2, length(rawdata), 3)]
  
  raw.counts$HybRefShort <- sub('\\|.*', '', rawdata$`Hybridization REF`)
  
  selected <- genes %>% inner_join(raw.counts, by = c('Gene.ID' = 'HybRefShort'))
  selected$HybRefShort <- NULL
  selected$Ensembl.ID <- NULL
  selected$`Hybridization REF` <- NULL
  
  rownames(selected) <- selected$Gene.ID
  
  selected$Gene.ID <- NULL
  
  counts <- as.data.frame(sapply(selected, as.integer))
  rownames(counts) <- rownames(selected)
  
  condition <- factor(barcode.to.type(colnames(counts)))
  
  counts <- counts[,condition != 'other']
  conditions <- factor(barcode.to.type(colnames(counts)))
  cds <- newCountDataSet(counts, conditions)
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  #plotDispEsts(cds)
  
  stats <- counts
  colnames(stats) <- unlist(barcode.to.type(colnames(stats)))
  stats.healthy <- stats[,grep('healthy', colnames(stats))]
  stats.tumor <- stats[,grep('tumor', colnames(stats))]
  
  stats$medianHealthy <- apply(stats.healthy, 1, median)
  stats$medianTumor <- apply(stats.tumor, 1, median)
  
  stats$IQRhealthy <- apply(stats.healthy, 1, IQR)
  stats$IQRtumor <- apply(stats.tumor, 1, IQR)
  
  stats$SEhealthy <- apply(stats.healthy, 1, se)
  stats$SEtumor <- apply(stats.tumor, 1, se)
  
  stats$Gene.ID <- rownames(stats)
  
  stats <- subset(stats, select = c('Gene.ID', 'medianHealthy', 'medianTumor', 'IQRhealthy', 'IQRtumor', 'SEhealthy', 'SEtumor'))
  
  res <- nbinomTest(cds, 'healthy', 'tumor')
  
  res <- res %>% inner_join(stats, by = c('id' = 'Gene.ID'))
  res <- rename(res, c('baseMeanA' = 'baseMeanHealthy', 'baseMeanB' = 'baseMeanTumor'))
  
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
