library(readr)
library(reshape2)

source('./data-utils.R')
source('./utils.R')
source('./subtypes/utils.R')

ignored <- c('age', 'OS_IND', 'OS_days', 'days_to_death', 'vital_status', 'days_to_last_followup', 
             'RFS_days', 'RFS_IND', 'TIME_TO_EVENT', 'days_to_new_tumor', 'number_pack_years_smoked', 
             'OS_months', 'RFS_months', 'Age', 'Days.to.Last.Followup', 'Disease.Free.Status', 
             'overall_survival_months', 'disease_free_months', 'pack_per_year',
             'days_to_last_follow_up', 'overall_survival_days')


mrna.dir <- '../../download/mrna-rsem-normalized-2015-11-01/'

files <- select.tumor.files(patient.freqs.path = '../../results/stats/normalized/stats-rsem-normalized-1-vs-11.csv',
                            files.path = '../../download/mrna-rsem-normalized-2015-11-01/')

subtypes.dir <- '../../KRAB_and_TCGA_Subtypes/Data v3 2016.06.06/'

res.path <- '../../../boxplots-subtypes-krabmen/resources'

tumor.names <- c()

for (subtype.file in list.files(subtypes.dir, pattern = 'csv', full.names = F)) {
  subtype.full.file <- file.path(subtypes.dir, subtype.file)
  
  subtypes.data <- read_delim(subtype.full.file, delim = ';', na = c('NA', 'na', ''))
  
  # Remove ignored
  subtypes.data <- subtypes.data[, setdiff(colnames(subtypes.data), ignored)]
  
  # Remove subtypes with only one value
  only.one <- apply(subtype.data, 2, function(x) length(unique(x))) == 1
  to.remove <- setdiff(names(only.one[only.one == T]), 'tumor')
  
  subtypes.data <- subtypes.data[, setdiff(colnames(subtypes.data), to.remove)]
  
  subtypes.data[is.na(subtypes.data)] <- 'N/A'
  
  tumor.name <- gsub('(.*)\\.csv', '\\1', subtype.file)
  
  # Find the file with tumor expression data
  tumor.file <- list.files(mrna.dir, pattern = tumor.name, recursive = T, full.names = T)
  
  results.file <- file.path(res.path, paste0(tumor.name, '.RData'))
  
  if (length(tumor.file) > 0 && !file.exists(results.file)) {
  
    print(tumor.name)
    
    file <- tumor.file
    
    tumor.data <- .load.tumor(tumor.name, file)
    
    join <- tumor.data %>% inner_join(subtypes.data, by = c('barcode' = 'sampleID'))
    
    join$barcode <- NULL
    
    long.format <- melt(join, 
                        id.vars = setdiff(colnames(subtypes.data), 'sampleID'), 
                        value.name = 'expression', 
                        variable.name = 'gene')
    
    long.format <- Filter(function(x) {
      (length(unique(x)) > 1)
    }, long.format)
    
    long.format$tumor <- tumor.name
    
    subtype.data <- long.format
    
    save(subtype.data, file = results.file)
    
    tumor.names <- c(tumor.names, tumor.name)
  }
}

write.table(tumor.names, 
            file = file.path(res.path, 'tumor.names.csv'), 
            row.names = F,
            col.names = c('name'))



# P-values
pval.results <- file.path(subtypes.dir, 'results')

for (tumor.results in list.files(pval.results, full.names = T)) {
  
  subtype.results <- list.files(tumor.results, full.names = T)
  
  tumor.name <- basename(tumor.results)

  print(paste('pvalues for ', tumor.name))
  
  if (tumor.name != 'README.md') {
    subt.tum.dir <- file.path(res.path, tumor.name)
    dir.create(subt.tum.dir, showWarnings = F)
    
    for (subtype.res in subtype.results) {
      raw <- read.csv(subtype.res)
      
      subtype <- gsub('(.*)\\.csv', '\\1', basename(subtype.res))
      
      # gene and p.vals, ignore mean, median, Q1, Q3, etc.
      pvals <- raw[, c(1, grep('p.adj', colnames(raw)))]
      
      save(pvals, 
           file = file.path(subt.tum.dir, paste0(subtype, '.RData')))
    }
  }
}
