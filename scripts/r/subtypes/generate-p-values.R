source('subtypes/utils.R')

subtypes.dir <- '../../KRAB_and_TCGA_Subtypes/Data v3 2016.06.06/'

mrna.dir <- '../../download/mrna-rsem-normalized-2015-11-01/'

# Ignored column names from KRAB_and_TCGA_Subtypes. 
# We don't perform analysis of variance for these subtypes.
ignored <- c('age', 'OS_IND', 'OS_days', 'days_to_death', 'vital_status', 'days_to_last_followup', 
             'RFS_days', 'RFS_IND', 'TIME_TO_EVENT', 'days_to_new_tumor', 'number_pack_years_smoked', 
             'OS_months', 'RFS_months', 'Age', 'Days.to.Last.Followup', 'Disease.Free.Status', 
             'overall_survival_months', 'disease_free_months', 'pack_per_year',
             'days_to_last_follow_up', 'overall_survival_days')

# Run to generate results, for each subtype file generate results
for (subtype.file in list.files(subtypes.dir, pattern = 'csv', full.names = F)) {
  
  #Extract name of the tumor
  tumor <- gsub('([a-zA-Z]+).*', '\\1', subtype.file)
  
  #Find the paht of the file with subtypes
  subtype.full.file <- file.path(subtypes.dir, subtype.file)
  
  # Find the file with tumor expression data
  tumor.file <- list.files(mrna.dir, pattern = tumor, recursive = T, full.names = T)
  
  # If file was found, perform analysis
  if (length(tumor.file) > 0) compare.subtypes(tumor, tumor.file, subtype.full.file, ignored)
}

