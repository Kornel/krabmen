library(ggplot2)
library(readr)
library(ggthemes)
library(dplyr)
library(reshape)

source('utils.R')

results.dir <- '../../results/heatmap'

genes <- read.csv('../../data/KRAB ZNF gene master list.csv', sep = '\t')
genes$Gene.ID <- as.character(genes$Gene.ID)

patient.freqs.path <- '../../results/stats/stats-rsem-normalized-1-vs-11.csv'
patient.freqs <- read.csv(patient.freqs.path) %>% filter(normal >= 9)

files.path <- '../../download/mrna-rsem-normalized-2015-11-01'

files <- list.files(path = files.path, pattern = '*-clean.txt', recursive = T)

filtered.files <- c()
  
for (tumor in patient.freqs$tumor.name) {
  contains <- lapply(files, function(x) grepl(tumor, x))
  file <- files[contains == T]
  filtered.files <- c(file, filtered.files)
}

missing.genes <- data.frame(tumor = character(0), gene.id = character(0))

heatmap.data <- data.frame(log2median = numeric(0), 
                           log2mean = numeric(0), 
                           Gene.Id = character(0), 
                           tumor.name = character(0))

for (file in filtered.files) {
  
  tumor.name <- sub('gdac.broadinstitute.org_(\\w*)\\..*', '\\1', file)
  
  print(sprintf("Reading %s", file))

  data <- read_delim(paste0(files.path, '/', file), delim = '\t')
  
  data$HybRefShort <- sub('\\|.*', '', data$`Hybridization REF`)
  
  joined <- genes %>% left_join(data, by = c('Gene.ID' = 'HybRefShort'))
    
  missing <- (joined %>% filter(is.na(`Hybridization REF`)))
  missing$tumor <- tumor.name

  missing.genes <- rbind(missing.genes, missing[, c('Gene.ID', 'tumor')])

  filtered <- joined %>% filter(!is.na(`Hybridization REF`))
  filtered$Ensembl.ID <- NULL
  filtered$`Hybridization REF` <- NULL

  t <- as.data.frame(filtered)
  rownames(t) <- t$Gene.ID
  t$Gene.ID <- NULL
  colnames(t) <- patient.code.to.type(colnames(t))
  
  t.tumor <- t[ , grep("tumor", colnames(t))]
  t.healthy <- t[ , grep("healthy", colnames(t))]
  
  t.tumor$median = apply(t.tumor, 1, safe.median)
  t.tumor$mean = apply(t.tumor, 1, safe.mean)
  
  t.healthy$median = apply(t.healthy, 1, safe.median)
  t.healthy$mean = apply(t.healthy, 1, safe.mean)
  
  t.tumor <- t.tumor[,c('median', 'mean')]
  t.healthy <- t.healthy[,c('median', 'mean')]
  
  colnames(t.tumor) <- paste0('tumor.', colnames(t.tumor))
  colnames(t.healthy) <- paste0('healthy.', colnames(t.healthy))
  
  t.tumor$Gene.Id <- rownames(t.tumor)
  t.healthy$Gene.Id <- rownames(t.healthy)
  
  t.joined <- t.tumor %>% inner_join(t.healthy, by = c('Gene.Id' = 'Gene.Id'))
  t.joined$log2median <- log2(t.joined$tumor.median/t.joined$healthy.median)
  t.joined$log2mean <- log2(t.joined$tumor.mean/t.joined$healthy.mean)
  
  t.joined$tumor.name <- tumor.name

  
  
  heatmap.data <- rbind(heatmap.data, t.joined)
}

write.csv(heatmap.data, file = sprintf('%s/heatmap-data.csv', results.dir), row.names = F)

library(pheatmap)

h <- heatmap.data[, c('log2median', 'log2mean', 'Gene.Id', 'tumor.name')]  


library(reshape2)
h.mean <- dcast(h, tumor.name ~ Gene.Id, value.var = 'log2mean')
rownames(h.mean) <- h.mean$tumor.name
h.mean$tumor.name <- NULL
h.mean <- t(h.mean)
pheatmap(h.mean, 
         show_colnames = T, 
         cluster_cols = T, 
         cluster_rows = T, 
         fontsize_row = 2,
         filename = sprintf('%s/heatmap-mean.png', results.dir),
         width = 10,
         height = 20)

h.median <- dcast(h, tumor.name ~ Gene.Id, value.var = 'log2median')
rownames(h.median) <- h.median$tumor.name
h.median$tumor.name <- NULL
h.median <- t(h.median)
pheatmap(h.median, 
         show_colnames = T, 
         cluster_cols = T, 
         cluster_rows = T, 
         fontsize_row = 2,
         filename = sprintf('%s/heatmap-median.png', results.dir),
         width = 10,
         height = 20)


# write.csv(missing.genes, file = sprintf('%s/missing-genes.csv', results.dir), row.names = F)
