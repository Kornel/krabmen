library(stringr)
library(plyr)

headers.dir <- '../../data/'
  
files <- list.files(headers.dir, recursive = T)
  
result <- data.frame(tumor.name = character(0), 
                     tumor = numeric(0), 
                     normal = numeric(0), 
                     stringsAsFactors = F)
  
results.dir <- '../../results'

for (f in files) {
  df <- read.csv(sprintf('%s/%s', headers.dir, f), nrows = 1, sep = '\t', header = F, stringsAsFactors = F)
  df$V1 <- NULL
    
  tumor.name <- sub('(\\w*)\\..*', '\\1', f)
    
  tumor <- 0
  normal <- 0
  
  participants <- c()
    
  for (x in df[1,]) {
    if (str_count(x, '-') != 6) {
      print(sprintf('Invalid barcode %s in file %s', x, tumor.name))  
    } else {
      elems <- unlist(strsplit(x, '-'))
      sample.vial <- elems[4]
      sample <- as.numeric(substr(sample.vial, 1, 2))
        
      participant <- elems[3]
        
      participants <- c(participant, participants)
        
      # https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
      # Tumor types range from 01 - 09, 
      # normal types from 10 - 19 and 
      # control samples from 20 - 29. 
        
      if (sample == 1) {
        tumor <- tumor + 1
      } else if (sample == 11) {
        normal <- normal + 1
      } 
        
    }
  }
  result[nrow(result) + 1, ] <- c(tumor.name, tumor, normal)

  dest.file <- sprintf('%s/table-%s.csv', results.dir, tumor.name)
  write.table(count(participants), file = dest.file, row.names = F)
}
  
dest.file <- sprintf('%s/stats-rsem-1-vs-11.csv', results.dir)
write.table(result, file = dest.file, row.names = F)

