library(ggplot2)
library(readr)
library(ggthemes)
library(dplyr)
library(reshape)

source('utils.R')

results.dir <- '../../results/boxplot'

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
other.types <- data.frame(tumor = character(0), barcode = character(0))

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

  t <- as.data.frame(t(filtered))
  colnames(t) <- unlist(as.list(t[1, ]))
  t <- t[-1,]
  t <- factorsNumeric(t)
  
  other <- c()
  
  t$type <- factor(unlist(lapply(rownames(t), patient.code.to.type)))
  
  t$barcodes <- rownames(t)
  
  other <- (t %>% filter(type == 'other'))$barcodes
  other.df <- data.frame(barcode = other)
  
  if (nrow(other.df) > 0) {
    other.df$tumor <- tumor.name
  }
  
  other.types <- rbind(other.types, other.df)
  
  t <- t %>% filter(type != 'other')
  
#   mdata <- melt(t, id = 'type')
#   
#   mdata <- mdata %>% filter(variable %in% c('ZNF195','TRIM28'))
  
  p <- ggplot(t, aes(y = ZNF195, x = type, fill = type)) + 
    geom_boxplot(alpha = 1, coef = 100 ) +
    scale_fill_manual(values = tableau_color_pal()(3)) + 
    ggtitle(tumor.name) #+ facet_wrap(~ variable, nrow = 1)
  
  ggsave(sprintf('%s/%s.png', results.dir, tumor.name), p, dpi = 180)
}

write.csv(other.types, file = sprintf('%s/other-types.csv', results.dir), row.names = F)
# 
# write.csv(missing.genes, file = sprintf('%s/missing-genes.csv', results.dir), row.names = F)
