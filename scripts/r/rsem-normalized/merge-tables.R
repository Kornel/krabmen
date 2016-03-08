library(readr)
library(data.table)

files <- list.files(path = '../../results/rsem-normalized/partial/', pattern = '*.csv', full.names = T)

all <- data.frame()

for (file in files) {
  tumor.name <- sub('.*/(\\w*)\\.csv', '\\1', file)
  
  print(tumor.name)
  
  raw <- read.csv(file)
  raw$tumor <- tumor.name
  
  all <- rbind(all, raw)
}

values <- setdiff(colnames(all), c('Gene.ID', 'tumor'))

all.wide <- as.data.frame(data.table::dcast(setDT(all), Gene.ID ~ tumor, value.var = values))

cnames <- colnames(all.wide)[-1]
order <- c(1, order(gsub(".*_(.*)", "\\1", cnames)) + 1)
cnames <- c('Gene.ID', cnames[order(gsub(".*_(.*)", "\\1", cnames))])

sorted <- all.wide[,cnames]

write.csv(all, file = '../../results/rsem-normalized/full-table-pvalues-long.csv', row.names = F)
write.csv(sorted, file = '../../results/rsem-normalized/full-table-pvalues.csv', row.names = F)