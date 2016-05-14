library(dplyr)
library(readr)
source('data-utils.R')
library(reshape2)

.load.tumor <- function(path) {
  
  print(paste('Loading tumor data from', path))
  
  prad <- read_delim(path, delim = '\t')
  
  raw.counts <- prad[-1,]
  raw.counts$HybRefShort <- hyb.ref.short(raw.counts)
  
  tumor.name <- 'PRAD'
  
  data <- get.gene.table(raw.counts, tumor.name, resolve.types = F)$gene.table 
  data$`Hybridization REF` <- NULL
  r <- rownames(data)
  data <- as.data.frame(sapply(data[,-1], as.numeric))
  rownames(data) <- r
  data <- as.data.frame(t(data))
  data$barcode <- substr(rownames(data), 1, 15)
  data <- data %>% filter(substr(barcode, 14, 15) == '01')
  data$barcode <- substr(data$barcode, 1, 12)
  data
}

.compare.subtypes <- function(subtypes.file, data) {
  
  genes <- setdiff(colnames(data), 'barcode')
  
  subtypes.file.name <- gsub('.*/(\\w+).csv', '\\1', subtypes.file)
  
  subtypes.data <- read_delim(subtypes.file, delim = ';', na = c('NA', 'na', ''))
  
  join <- data %>% inner_join(subtypes.data, by = c('barcode' = 'sampleID'))
  
  subtypes <- setdiff(colnames(join), c(genes, 'barcode'))
  
  ignored <- c('age', 'OS_IND', 'OS_days', 'days_to_death', 'vital_status', 'days_to_last_followup', 
               'RFS_days', 'RFS_IND', 'TIME_TO_EVENT', 'days_to_new_tumor', 'number_pack_years_smoked', 
               'OS_months', 'RFS_months', 'Age', 'Days.to.Last.Followup')
  
  subtypes <- setdiff(subtypes, ignored)
  
  for (subtype in subtypes) {
    
    na <- is.na(join[,subtype])
    
    subtCol <- as.character(join[,subtype])[!na]
    
    subtypes.str <- toString(paste(unlist(unique(subtCol)), collapse = ', '))
    
    if (length(unique(subtCol)) <= 1) {
      print(sprintf('Cannot perform comparison, only values of %s found: %s', subtype, subtypes.str))
    } else {
      print(sprintf('For %s found: %s', subtype, subtypes.str))
      
      dir <- file.path('../../KRAB_and_TCGA_Subtypes/results', subtypes.file.name)
      dir.create(dir, showWarnings = F, recursive = T)
      file <- file.path(dir, paste0(subtype, '.csv'))
      
      b <- data.frame()
      for (gene in genes) {
        
        geneCol <- join[,gene][!na]
        a <- aov(geneCol ~ subtCol)
        t.hsd <- as.data.frame(TukeyHSD(a)$`subtCol`)
        t.hsd$type <- rownames(t.hsd)
        t.hsd$gene <- gene
        b <- rbind(b, dcast(t.hsd, gene ~ type, value.var = c('p adj')))
      }
      write.table(b, file = file, row.names = F, sep = ',')
    }
  }
}

compare.subtypes <- function(tumor.file, subtypes.file) {
  data <- .load.tumor(tumor.file)
  .compare.subtypes(subtypes.file, data)
}