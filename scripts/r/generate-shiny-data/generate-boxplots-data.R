library(readr)
library(reshape2)

source('./data-utils.R')
source('./utils.R')

files <- select.tumor.files(patient.freqs.path = '../../results/stats/normalized/stats-rsem-normalized-1-vs-11.csv',
                            files.path = '../../download/mrna-rsem-normalized-2015-11-01/')

file <- files[1]

all.tumors <- data.frame()
  
for (file in files) {
  tumor.name <- filename.to.tumor(file)
  
  print(tumor.name)
  rawdata <- read_delim(file, delim = '\t')  
  rawdata$HybRefShort <- hyb.ref.short(rawdata)
  rawdata$`Hybridization REF` <- NULL
  gene.table <- get.gene.table(rawdata, tumor.name)$gene.table
  
  gene.table <- t(gene.table)
  r <- rownames(gene.table)
  rownames(gene.table) <- NULL
  gene.table <- as.data.frame(gene.table)
  gene.table$type <- r
  
  gene.table <- gene.table %>% filter(type != 'other')
  
  gene.table$tumor <- tumor.name

  all.tumors <- rbind(all.tumors, gene.table)
  # save(gene.table, file = paste0('../../../boxplots-krabmen/res2/', tumor.name, '.Rdata'))  
}

# {
#   lung <- read_delim('../../download/mrna-rsem-normalized-2015-11-01/LUXX/LUXX.csv', delim = '\t')  
#   lung$HybRefShort <- hyb.ref.short(lung)
#   lung$`Hybridization REF` <- NULL
#   gene.table <- get.gene.table(lung, tumor.name)$gene.table
#   
#   gene.table <- t(gene.table)
#   r <- rownames(gene.table)
#   rownames(gene.table) <- NULL
#   gene.table <- as.data.frame(gene.table)
#   gene.table$type <- r
#   
#   gene.table <- gene.table %>% filter(type != 'other')
#   
#   gene.table$tumor <- 'LUNG'
#   
#   all.tumors <- rbind(all.tumors, gene.table)
# }
# 




all.tumors.long <- melt(all.tumors, id.vars = c('tumor', 'type'), value.name = 'expression', variable.name = 'gene')
save(all.tumors.long, file = paste0('../../../boxplots-krabmen/resources/all.Rdata'))

write.table(unique(all.tumors$tumor), file = paste0('../../../boxplots-krabmen/resources/tumor-names.csv'), row.names = F,
          col.names = c('name'), sep = ',')

write.table(setdiff(colnames(all.tumors), c('tumor', 'type')), file = paste0('../../../boxplots-krabmen/resources/gene-names.csv'), row.names = F,
            col.names = c('name'), sep = ',')


pvalues <- read_csv('../../results/bionmialTest/full-table-pvalues-long.csv')
pvalues <- subset(pvalues, select = c('pval', 'padj', 'id', 'tumor'))
colnames(pvalues)[3] <- 'gene'
save(pvalues, file = '../../../boxplots-krabmen/resources/pvalues.Rdata')
