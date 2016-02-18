library(readr)

files.path <- '../../download/mrna-rsem-normalized-2015-11-01'

files <- list.files(path = files.path, pattern = '*data.txt', recursive = T)
files2 <- head(files, 1)

dest.dir <- '../../results/datasets/cancer-only/'

for (file in files) {
  
  tumor.name <- sub('.*?/(\\w*)\\..*', '\\1', file)
  
  print(tumor.name)
  
  df <- read_delim(sprintf('%s/%s', files.path, file), delim = '\t')
  df <- df[-c(1), ] # First row contains comments (gene_id, normalized_count, ...)
  
  mask <- sapply(colnames(df)[-1], USE.NAMES = F   ,  function(x) {
    elems <- unlist(strsplit(x, '-'))
    sample.vial <- elems[4]
    sample <- as.numeric(substr(sample.vial, 1, 2))
    if (sample == 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  
  mask <- c(TRUE, mask) #Take first column (gene names) as well

  df <- df[, mask]
  
  write.csv(df, file = sprintf('%s/%s.csv', dest.dir, tumor.name), row.names = F)
  
}
