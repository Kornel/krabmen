library(DESeq)
library(dplyr)
library(readr)

source('./data-utils.R')
source('./utils.R')
source('./binomtest/stats.R')

lung.path <- '../../download/mrna-rsem-2015-11-01/LUNG/LUNG.csv'

# Merge LUAD and LUSC into LUNG
# if (!file.exists(luxx.path)) {

if (TRUE) {
  luad.path <- '../../download/mrna-rsem-2015-11-01/gdac.broadinstitute.org_LUAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2015110100.0.0/LUAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt.gz'
  lusc.path <- '../../download/mrna-rsem-2015-11-01/gdac.broadinstitute.org_LUSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2015110100.0.0/LUSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt.gz'
  
  rawdata.luad <- read_delim(luad.path, delim = '\t')
  rawdata.lusc <- read_delim(lusc.path, delim = '\t')
  
  raw.counts.luad <- raw.counts.only(rawdata.luad)
  raw.counts.lusc <- raw.counts.only(rawdata.lusc)
  
  filtered.counts.luad <- get.gene.table(raw.counts.luad, 'LUAD')$gene.table
  filtered.counts.lusc <- get.gene.table(raw.counts.lusc, 'LUSC')$gene.table
  
  filtered.counts.luad$id = rownames(filtered.counts.luad)
  filtered.counts.lusc$id = rownames(filtered.counts.lusc)
  
  n <- 0
  L <- length(colnames(filtered.counts.luad)) - 1
  colnames(filtered.counts.luad)[1:L] <- sapply(colnames(filtered.counts.luad)[1:L],
         function(x) {
           n <<- n + 1
           return(paste(x, n, sep = '-'))
         }
  )
  
  n <- L
  L <- length(colnames(filtered.counts.lusc)) - 1
  colnames(filtered.counts.lusc)[1:L] <- sapply(colnames(filtered.counts.lusc)[1:L],
                                                function(x) {
                                                  n <<- n + 1
                                                  return(paste(x, n, sep = '-'))
                                                }
  )
  
  lung <- filtered.counts.luad %>% inner_join(filtered.counts.lusc, by=c('id' = 'id'))
  
  rownames(lung) <- lung$id
  lung$id <- NULL
  
  colnames(lung) <- sapply(colnames(lung), function(x) {
    return(unlist(strsplit(x, '-'))[[1]])
  })
  
  
  tumor.name <- 'LUNG'
  
  counts <- as.data.frame(sapply(lung, as.integer))
  rownames(counts) <- rownames(lung)
  
  counts <- counts[, colnames(counts) != 'other']
  conditions <- factor(gsub('(.*)\\..*', '\\1', colnames(counts)))
  
  cds <- newCountDataSet(counts, conditions)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  res <- nbinomTest(cds, 'healthy', 'tumor')
  
  res <- res %>% inner_join(get.stats(counts), by = c('id' = 'Gene.ID'))
  res <- plyr::rename(res, c('baseMeanA' = 'baseMeanHealthy', 'baseMeanB' = 'baseMeanTumor'))
  
  write.csv(res, file = sprintf('../../results/bionmialTest/partial/%s.csv', tumor.name), row.names = F)
}

