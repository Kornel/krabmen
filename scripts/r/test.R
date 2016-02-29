library(readr)
library(edgeR)
library(dplyr)

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

tissue <- barcode.to.type(colnames(rawdata)[-1])

y <- DGEList(counts = rawdata[, 2:ncol(rawdata)], genes = rawdata[,1], group = tissue)
y1 <- calcNormFactors(y)

# y1$sample
# 
samples <- colnames(y1)
tissue <- barcode.to.type(samples)
patient <- barcode.to.participant(samples)

design <- model.matrix(~patient + tissue)

y2 <- estimateDisp(y1, design, robust = F)

fit <- glmFit(y1, design)
