library(dplyr)

get.gene.table <- function(data, tumor.name) {
 
  genes <- read.csv('../../data/KRAB ZNF gene master list.csv', sep = '\t')
  genes$Gene.ID <- as.character(genes$Gene.ID)
  
  joined <- genes %>% inner_join(data, by = c('Gene.ID' = 'HybRefShort'))
  missing <- data.frame(Gene.ID = setdiff(genes$Gene.ID, data$HybRefShort), tumor = tumor.name)
  joined$Ensembl.ID <- NULL  
  
  rownames(joined) <- joined$Gene.ID
  joined$Gene.ID <- NULL
  colnames(joined) <- patient.code.to.type(colnames(joined))
  
  return(list(gene.table = as.data.frame(joined), missing = missing))
}

raw.counts.only <- function(rawdata) {
  data <- rawdata[,seq(2, length(rawdata), 3)]
  data$HybRefShort <- sub('\\|.*', '', rawdata$`Hybridization REF`)
  return(data)
}

scaled.only <- function(rawdata) {
  data <- rawdata[,seq(3, length(rawdata), 3)]
  data$HybRefShort <- sub('\\|.*', '', rawdata$`Hybridization REF`)
  return(data)
}


select.tumor.files <- function(patient.freqs.path = '../../results/stats/raw/stats-rsem-1-vs-11.csv',
                               files.path = '../../download/mrna-rsem-2015-11-01/') {
  
  results.dir <- '../../results/heatmap'
  
  patient.freqs.path <- '../../results/stats/raw/stats-rsem-1-vs-11.csv'
  patient.freqs <- read.csv(patient.freqs.path) %>% filter(normal >= 9)
  
  files <- list.files(path = files.path, pattern = '*-clean.txt', recursive = T, full.names = T)
  
  filtered.files <- c()
  
  for (tumor in patient.freqs$tumor.name) {
    contains <- lapply(files, function(x) grepl(tumor, x))
    file <- files[contains == T]
    filtered.files <- c(file, filtered.files)
  }
  
  return(filtered.files)
}