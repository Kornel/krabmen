library(DESeq)
library(dplyr)
library(readr)

genes <- read.csv('../../data/KRAB ZNF gene master list.csv', sep = '\t')
genes$Gene.ID <- as.character(genes$Gene.ID)

source('utils.R')

path <- '../../download/mrna-rsem-normalized-2015-11-01/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt-clean.txt'

rawdata <- read_delim(file = path, delim = '\t')
rawdata$HybRefShort <- sub('\\|.*', '', rawdata$`Hybridization REF`)

rawdata <- genes %>% inner_join(rawdata, by = c('Gene.ID' = 'HybRefShort'))
rawdata$HybRefShort <- NULL
rawdata$Ensembl.ID <- NULL
rawdata$`Hybridization REF` <- NULL

countTable <- rawdata[,-1]
rownames(countTable) <- rawdata$Gene.ID

condition <- factor(barcode.to.type(colnames(rawdata)[-1]))

count2 <- countTable[,condition != 'other']
cond2 <- factor(barcode.to.type(colnames(count2)))
cds <- newCountDataSet(count2, cond2)

