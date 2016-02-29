# if (!require(devtools)) {
#   install.packages("devtools")
#   require(devtools)
# }
# 
# install_github("RTCGA/RTCGA")

library(RTCGA)

cohorts <- infoTCGA() %>% rownames() %>% sub('-counts', '', x=.)

release.dates <- c('2015-11-01')

download.dir <- '../../download/'


  for (release.date in release.dates) {
    dest.dir <- sprintf('%s/mrna-rsem-%s', download.dir, release.date)
    
    dir.create(dest.dir, showWarnings = F)
    
    sapply(cohorts, function(element){
      tryCatch({
        downloadTCGA(cancerTypes = element, 
                     dataSet ='RSEM_genes__data.Level_3',
                     destDir = dest.dir, 
                     date = release.date)
      },
      error = function(cond){
        cat("Error: Maybe there weren't mutations data for ", element, " cancer.\n")
      }
      )
    })
    
  }