library(readr)
library(dplyr)
library(IlluminaHumanMethylation27k.db)
library(reshape2)
source('subtypes/utils.R')
brca.file <- '../../download/mrna-rsem-normalized-2015-11-01/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015110100.0.0/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt.gz'
subtype.file <- 'methylation/methyl-clusters.csv'
compare.subtypes('BRCA', brca.file, subtype.file, c(''), only.krab = F)

brca <- read_delim('methylation/results/methylation/methyl-clusters.csv/methylation.Clusters.all.csv', delim=',')
mapping450k <- read.delim('../../data/IlluminaHumanMethylation450k-probes-mapping.csv.gz', sep = ',', col.names = c('probe', 'gene', 'multi'))

ProbeToSymbol <- IlluminaHumanMethylation27kSYMBOL
mapped_probes <- mappedkeys(ProbeToSymbol)
mapped_probes.df <- as.data.frame(ProbeToSymbol[mapped_probes])

brca$probe27k <- unname(sapply(brca$gene, function(gene) {
  symbols <- mapped_probes.df[mapped_probes.df$symbol == gene, 1]
  paste(sort(symbols), collapse=' ')
}))

mapping450k$gene <- as.character(mapping450k$gene)

s2eg <- org.Hs.egSYMBOL2EG
mapping <- as.data.frame(s2eg[mappedkeys(s2eg)])
probe2symbol <- mapping450k %>% 
  inner_join(mapping, by = c('gene' = 'gene_id')) %>%
  select(probe, symbol) %>%
  rename(probe450k = probe)

brca$probe450k <- unname(sapply(brca$gene, function(gene) {
  probes <- probe2symbol[probe2symbol$symbol == gene,1]
  paste(sort(probes), collapse=' ')
}))

write.table(brca, file = 'methylation/results/BRCA-probes-mapped.csv', row.names = F, sep = ',')
