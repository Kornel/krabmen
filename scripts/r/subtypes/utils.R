library(dplyr)
library(readr)
library(reshape2)

# Load tumor data
.load.tumor <- function(tumor.name, path) {
  
  print(paste('Loading tumor data from', path))
  
  raw <- read_delim(path, delim = '\t')
  
  # Ignore first column
  raw.counts <- raw[-1,]
  
  # Shorten gene names
  raw.counts$HybRefShort <- .hyb.ref.short(raw.counts)
  
  # Filter for relevant genes
  data <- .get.gene.table(raw.counts, tumor.name) 
  
  #Fetch gene names
  r <- rownames(data)
  
  # Convert to data to numbers (first column removed earlier caused the data frame to be loaded as non numeric data)
  data <- as.data.frame(sapply(data[,-1], as.numeric))
  rownames(data) <- r
  
  #Transpose, so we have genes in columns and patients in rows
  data <- as.data.frame(t(data))
  
  #Shorten barcode
  data$barcode <- substr(rownames(data), 1, 15)
  
  #Filter only relevant patients with tumor samples
  data <- data %>% filter(substr(barcode, 14, 15) == '01')
  data$barcode <- substr(data$barcode, 1, 12)
  data
}

.compare.subtypes <- function(tumor.name, subtypes.file, data, ignored) {
  
  genes <- setdiff(colnames(data), 'barcode')
  
  subtypes.file.name <- gsub('.*/(\\w+).csv', '\\1', subtypes.file)
  
  subtypes.data <- read_delim(subtypes.file, delim = ';', na = c('NA', 'na', ''))
  
  join <- data %>% inner_join(subtypes.data, by = c('barcode' = 'sampleID'))
  
  subtypes <- setdiff(colnames(join), c(genes, 'barcode'))
  
  subtypes <- setdiff(subtypes, ignored)
  
  # For each subtype perform analysis
  for (subtype in na.omit(subtypes)) {
  
    na <- is.na(join[,subtype])
    
    subtCol <- as.character(join[,subtype])[!na]
    
    subtypes.str <- toString(paste(unlist(unique(subtCol)), collapse = ', '))
    
    # If there is actually more than one value in subtype column, analyse
    if (length(unique(subtCol)) <= 1) {
      print(sprintf('IGNORED subtype \'%s\'. Cannot perform comparison, only values found: %s', subtype, subtypes.str))
    } else {
      print(sprintf('OK subtype \'%s\'. Values found: %s', subtype, subtypes.str))
      
      dir <- file.path(dirname(subtypes.file), 'results', subtypes.file.name)
      dir.create(dir, showWarnings = F, recursive = T)
      file <- file.path(dir, paste0(subtype, '.csv'))
      
      b <- data.frame()
      
      # Perform comparison for each gene
      for (gene in genes) {
        
        geneCol <- join[,gene][!na]
        
        a <- aov(geneCol ~ subtCol)
        t.hsd <- as.data.frame(TukeyHSD(a)$`subtCol`)
        
        t.hsd$type <- rownames(t.hsd)
        t.hsd$gene <- gene
        t.hds.df <- dcast(t.hsd, gene ~ type, value.var = c('p adj'))
        colnames(t.hds.df)[-1] <- paste0('p-adj-', colnames(t.hds.df)[-1])
        
        t.hds.df <- .append.stat(geneCol, subtCol, mean, 'mean-', t.hds.df)
        t.hds.df <- .append.stat(geneCol, subtCol, function(x) sqrt(var(x) / length(x)), 'SEM-', t.hds.df)
        t.hds.df <- .append.stat(geneCol, subtCol, function(x) quantile(x, 1/4), 'Q1-', t.hds.df)
        t.hds.df <- .append.stat(geneCol, subtCol, function(x) quantile(x, 3/4), 'Q3-', t.hds.df)
        t.hds.df <- .append.stat(geneCol, subtCol, length, 'count-', t.hds.df)
        
        b <- rbind(b, t.hds.df)
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

# Extract gene name
.hyb.ref.short <- function(rawdata) {
  sub('\\|.*', '', rawdata$`Hybridization REF`)
}

# Filter relevant genes
.get.gene.table <- function(rawdata, tumor.name) {
  
  genes <- read.csv('../../data/KRAB ZNF gene master list.csv', sep = '\t')
  genes$Gene.ID <- as.character(genes$Gene.ID)
  
  joined <- genes %>% inner_join(rawdata, by = c('Gene.ID' = 'HybRefShort'))
  
  joined$Ensembl.ID <- NULL  
  
  rownames(joined) <- joined$Gene.ID
  joined$Gene.ID <- NULL
  joined$`Hybridization REF` <- NULL
  
  as.data.frame(joined)
  
}

compare.subtypes <- function(tumor.name, tumor.file, subtypes.file, ignored) {
  data <- .load.tumor(tumor.name, tumor.file)
  .compare.subtypes(tumor.name, subtypes.file, data, ignored)
}