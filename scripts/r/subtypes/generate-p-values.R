source('subtypes/utils.R')

subtypes.dir <- '../../KRAB_and_TCGA_Subtypes/Data v3 2016.06.06/'

mrna.dir <- '../../download/mrna-rsem-normalized-2015-11-01/'

# Merge LUAD and LUSC into LUXX

luxx.path <- '../../download/mrna-rsem-normalized-2015-11-01/LUXX/LUXX.csv'

if (!file.exists(luxx.path)) {
  
  luad.path <- '../../download/mrna-rsem-normalized-2015-11-01/gdac.broadinstitute.org_LUAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/LUAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt.gz'
  lusc.path <- '../../download/mrna-rsem-normalized-2015-11-01/gdac.broadinstitute.org_LUSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/LUSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt.gz'
  
  
  luad <- read_delim(luad.path, delim = '\t')
  lusc <- read_delim(lusc.path, delim = '\t')
  luxx <- luad %>% inner_join(lusc, by=c('Hybridization REF' = 'Hybridization REF'))
  
  write_delim(luxx, path = luxx.path, delim = '\t')
}

# Merge subtypes

# Ignored column names from KRAB_and_TCGA_Subtypes. 
# We don't perform analysis of variance for these subtypes.
ignored <- c('age', 'OS_IND', 'OS_days', 'days_to_death', 'vital_status', 'days_to_last_followup', 
             'RFS_days', 'RFS_IND', 'TIME_TO_EVENT', 'days_to_new_tumor', 'number_pack_years_smoked', 
             'OS_months', 'RFS_months', 'Age', 'Days.to.Last.Followup', 'Disease.Free.Status', 
             'overall_survival_months', 'disease_free_months', 'pack_per_year',
             'days_to_last_follow_up', 'overall_survival_days')

# Run to generate results, for each subtype file generate results
# for (subtype.file in list.files(subtypes.dir, pattern = 'csv', full.names = F)) {
for (subtype.file in c('LUXX.csv')) {
  
  #Extract name of the tumor
  tumor <- gsub('([a-zA-Z]+).*', '\\1', subtype.file)
  
  #Find the paht of the file with subtypes
  subtype.full.file <- file.path(subtypes.dir, subtype.file)
  
  # Find the file with tumor expression data
  tumor.file <- list.files(mrna.dir, pattern = tumor, recursive = T, full.names = T)
  
  # If file was found, perform analysis
  if (length(tumor.file) > 0) compare.subtypes(tumor, tumor.file, subtype.full.file, ignored)
}

