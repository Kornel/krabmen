library(dplyr)
library(readr)
source('data-utils.R')
source('utils.R')
library(reshape2)

.load.tumor <- function(tumor.name, path) {
  
  print(paste('Loading tumor data from', path))
  
  raw <- read_delim(path, delim = '\t')
  
  raw.counts <- raw[-1,]
  raw.counts$HybRefShort <- hyb.ref.short(raw.counts)
  
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

.compare.subtypes <- function(tumor.stats, tumor.name, subtypes.file, data, ignored) {
  
  genes <- setdiff(colnames(data), 'barcode')
  
  subtypes.file.name <- gsub('.*/(\\w+).csv', '\\1', subtypes.file)
  
  subtypes.data <- read_delim(subtypes.file, delim = ';', na = c('NA', 'na', ''))
  
  join <- data %>% inner_join(subtypes.data, by = c('barcode' = 'sampleID'))
  
  subtypes <- setdiff(colnames(join), c(genes, 'barcode'))
  
  subtypes <- setdiff(subtypes, ignored)
  
  for (subtype in na.omit(subtypes)) {
  
    na <- is.na(join[,subtype])
    
    subtCol <- as.character(join[,subtype])[!na]
    
    subtypes.str <- toString(paste(unlist(unique(subtCol)), collapse = ', '))
    
    if (length(unique(subtCol)) <= 1) {
      print(sprintf('IGNORED subtype \'%s\'. Cannot perform comparison, only values found: %s', subtype, subtypes.str))
    } else {
      print(sprintf('OK subtype \'%s\'. Values found: %s', subtype, subtypes.str))
      
      dir <- file.path('../../KRAB_and_TCGA_Subtypes/results', subtypes.file.name)
      dir.create(dir, showWarnings = F, recursive = T)
      file <- file.path(dir, paste0(subtype, '.csv'))
      
      b <- data.frame()
      for (gene in genes) {
        
        extra.stats <- tumor.stats[as.character(full.stats$Gene.ID) == gene,]
        extra.stats$Gene.ID <- NULL
        
        geneCol <- join[,gene][!na]
        
        a <- aov(geneCol ~ subtCol)
        t.hsd <- as.data.frame(TukeyHSD(a)$`subtCol`)
        
        t.hsd$type <- rownames(t.hsd)
        t.hsd$gene <- gene
        t.hds.df <- dcast(t.hsd, gene ~ type, value.var = c('p adj'))
        colnames(t.hds.df) <- paste0('p-adj', colnames(t.hds.df))
        
        t.hds.df <- .append.stat(geneCol, subtCol, mean, 'mean-', t.hds.df)
        t.hds.df <- .append.stat(geneCol, subtCol, se, 'SEM-', t.hds.df)
        t.hds.df <- .append.stat(geneCol, subtCol, function(x) quantile(x, 1/4), 'Q1-', t.hds.df)
        t.hds.df <- .append.stat(geneCol, subtCol, function(x) quantile(x, 3/4), 'Q3-', t.hds.df)
        
        x <- cbind(t.hds.df, extra.stats)
        
        b <- rbind(b, x)
      }
      write.table(b, file = file, row.names = F, sep = ',')
      print(sprintf('Analysis written to %s', file))
    }
  }
}

.append.stat <- function(geneCol, subtCol, stat, name, df) {
  stats <- tapply(geneCol, subtCol, stat)
  names(stats) <- paste0(name, names(stats))
  cbind(df, t(as.data.frame(stats)))
}

compare.subtypes <- function(tumor.stats, tumor.name, tumor.file, subtypes.file, ignored) {
  data <- .load.tumor(tumor.name, tumor.file)
  .compare.subtypes(tumor.stats, tumor.name, subtypes.file, data, ignored)
}