library(stringr)
library(plyr)

headers.dir <- '../../data/headers/rsem-normalized/'
  
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
  
  counts <- list()
  
  for (x in df[1,]) {
    if (str_count(x, '-') != 6) {
      print(sprintf('Invalid barcode %s in file %s', x, tumor.name))  
    } else {
      elems <- unlist(strsplit(x, '-'))
      sample.vial <- elems[4]
      sample <- as.numeric(substr(sample.vial, 1, 2))
        
      participant <- elems[3]
      
      if (is.null(counts[[participant]]$Tumor)) counts[[participant]]$Tumor <- 0
      if (is.null(counts[[participant]]$Normal)) counts[[participant]]$Normal <- 0
        
      # https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
      # Tumor types range from 01 - 09, 
      # normal types from 10 - 19 and 
      # control samples from 20 - 29. 
        
      if (sample == 1) {
        tumor <- tumor + 1
        counts[[participant]]$Tumor <- counts[[participant]]$Tumor + 1
      } else if (sample == 11) {
        normal <- normal + 1
        counts[[participant]]$Normal <- counts[[participant]]$Normal + 1
      } 
      
    }
  }
  result[nrow(result) + 1, ] <- c(tumor.name, tumor, normal)

  dest.file <- sprintf('%s/table-%s.csv', results.dir, tumor.name)
  counts.df <- do.call("rbind", lapply(counts, as.data.frame))
  counts.df <- cbind(Patient=rownames(counts.df), counts.df)
  write.table(counts.df, file = dest.file, row.names = F, sep = ",")
}
  
dest.file <- sprintf('%s/stats-rsem-1-vs-11.csv', results.dir)
write.table(result, file = dest.file, row.names = F, sep = ",")

