# if (!require(devtools)) {
#   install.packages("devtools")
#   require(devtools)
# }
# 
# install_github("RTCGA/RTCGA")

library(RTCGA)

cohorts <- infoTCGA() %>% rownames() %>% sub('-counts', '', x=.)

# cohorts <- c('STES', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM')
cohorts <- c('COAD','COADREAD','READ','UCEC')

release.dates <- c('2015-11-01')

download.dir <- '../../download/'

for (release.date in release.dates) {
  dest.dir <- sprintf('%s/mrna-rsem-%s', download.dir, release.date)
  dest.dir.norm <- sprintf('%s/mrna-rsem-normalized-%s', download.dir, release.date)
    
  dir.create(dest.dir, showWarnings = F)
    
  sapply(cohorts, function(element){
    tryCatch({
      downloadTCGA(cancerTypes = element, 
                   dataSet ='Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3',
                   destDir = dest.dir, 
                   date = release.date)
      
#       downloadTCGA(cancerTypes = element, 
#                    dataSet ='Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data',
#                    destDir = dest.dir.norm, 
#                    date = release.date)
    },
    error = function(cond){
      print(cond)
    }
    )
  })
}
