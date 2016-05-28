source('subtypes/utils.R')

subtypes.dir <- '../../KRAB_and_TCGA_Subtypes/'
mrna.dir <- '../../download/mrna-rsem-normalized-2015-11-01/'

# Ignored column names from KRAB_and_TCGA_Subtypes
ignored <- c('age', 'OS_IND', 'OS_days', 'days_to_death', 'vital_status', 'days_to_last_followup', 
             'RFS_days', 'RFS_IND', 'TIME_TO_EVENT', 'days_to_new_tumor', 'number_pack_years_smoked', 
             'OS_months', 'RFS_months', 'Age', 'Days.to.Last.Followup')

# Generated using rsem-normalized/generate-tables.R and rsem-normalized/merge-tables.R
full.stats <- read.csv('../../results/rsem-normalized/full-table-long.csv')

# Run to generate results
for (subtype.file in list.files(subtypes.dir, pattern = 'csv', full.names = F)) {
  tumor <- gsub('([a-zA-Z]+).*', '\\1', subtype.file)
  
  subtype.full.file <- file.path(subtypes.dir, subtype.file)
  
  tumor.file <- list.files(mrna.dir, pattern = tumor, recursive = T, full.names = T)
  
  tumor.stats <- full.stats[as.character(full.stats$tumor) == tumor,]
  
  if (length(tumor.file) > 0) compare.subtypes(tumor.stats, tumor, tumor.file, subtype.full.file, ignored)
}

