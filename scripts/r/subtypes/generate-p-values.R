source('subtypes/utils.R')

subtypes.dir <- '../../KRAB_and_TCGA_Subtypes/'
mrna.dir <- '../../download/mrna-rsem-normalized-2015-11-01/'

for (subtype.file in list.files(subtypes.dir, pattern = 'csv', full.names = F)) {
  tumor <- gsub('([a-zA-Z]+).*', '\\1', subtype.file)
  subtype.full.file <- file.path(subtypes.dir, subtype.file)
  
  tumor.file <- list.files(mrna.dir, pattern = tumor, recursive = T, full.names = T)
  
  if (length(tumor.file) > 0) compare.subtypes(tumor.file, subtype.full.file)
}
