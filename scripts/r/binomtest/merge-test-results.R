library(dplyr)
library(readr)
library(reshape2)
library(data.table)

files <- list.files(path = '../../results/bionmialTest/partial/',
                    pattern = '*.csv',
                    full.names = T)

all <- data.frame()

# file <- '../../results/bionmialTest/partial/LUSC.csv'

for (file in files) {

  tumor.name <- sub('.*/(\\w*)\\.csv', '\\1', file)

  print(tumor.name)

  raw <- read.csv(file, stringsAsFactors = F)
  raw$tumor <- tumor.name

  all <- rbind(all, raw)
}

values <- setdiff(colnames(all), c('id', 'tumor'))

all.wide <- as.data.frame(data.table::dcast(setDT(all), id ~ tumor, value.var = values))

cnames <- colnames(all.wide)[-1]
order <- c(1, order(gsub(".*_(.*)", "\\1", cnames)) + 1)

cnames <- c('id', cnames[order(gsub(".*_(.*)", "\\1", cnames))])

sorted <- all.wide[,cnames]

write.csv(all, file = '../../results/bionmialTest/full-table-pvalues-long.csv', row.names = F)
write.csv(sorted, file = '../../results/bionmialTest/full-table-pvalues.csv', row.names = F)
