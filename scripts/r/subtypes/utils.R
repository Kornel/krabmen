library(dplyr)
source('data-utils.R')

load.tumor <- function(path) {
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
  data
}

generate.pvalues <- function(data, subtypes.file, output) {
  
  subtype.mapping <- read.csv('../../data/LUSC.clinical.basic.tsv', sep = '\t', stringsAsFactors = F, col.names = c('barcode', 'subtype')) %>% filter(nchar(subtype) > 0)
  subtype.mapping$subtype <- as.character(subtype.mapping$subtype)
  gene.list <- colnames(data)[-length(colnames(data))]
  
  subtype.data <- data %>% inner_join(subtype.mapping, by = c('barcode' = 'barcode'))
  
  if (length(unique(subtype.data$subtype)) == 1) {
    print(sprintf('Cannot perform comparison, only one subtype found %s', toString(unique(subtype.data$subtype))))
  } else {
    b <- data.frame()
    for (gene in gene.list) {
      a <- aov(subtype.data[,gene] ~ subtype.data$subtype)
      t.hsd <- as.data.frame(TukeyHSD(a)$`subtype.data$subtype`)
      t.hsd$type <- rownames(t.hsd)
      t.hsd$gene <- gene
      b <- rbind(b, dcast(t.hsd, gene ~ type, value.var = c('p adj')))
    }
    
    write.table(b, file = output, row.names = F, sep = ',')
  }
}