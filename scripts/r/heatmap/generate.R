library(ggplot2)
library(readr)
library(ggthemes)
library(dplyr)
library(reshape)

source('utils.R')
source('data-utils.R')
source('heatmap-prepare.R')

missing.genes <- data.frame(tumor = character(0), gene.id = character(0))

heatmap.data <- data.frame(log2median = numeric(0), 
                           log2mean = numeric(0), 
                           Gene.Id = character(0), 
                           tumor.name = character(0))

files <- select.tumor.files(patient.freqs.path = '../../results/stats/normalized/stats-rsem-normalized-1-vs-11.csv',
                            files.path = '../../download/mrna-rsem-normalized-2015-11-01/')

for (file in files) {
  
  tumor.name <- sub('.*gdac.broadinstitute.org_(\\w*)\\..*', '\\1', file)
  
  print(sprintf("Reading %s", tumor.name))

  rawdata <- read_delim(file, delim = '\t')
  data <- rawdata
  data$`Hybridization REF` <- NULL
  data$HybRefShort <- sub('\\|.*', '', rawdata$`Hybridization REF`)
  
  gene.table.result <- get.gene.table(data, tumor.name)
  gene.table <- gene.table.result$gene.table
  
  missing.genes <- rbind(missing.genes, gene.table.result$missing)
  
  h <- prepare.heatmap.data(gene.table)  
  heatmap.data <- rbind(heatmap.data, h)
}

results.dir <- '../../results/heatmap/'

write.csv(heatmap.data, file = sprintf('%s/heatmap-data.csv', results.dir), row.names = F)
