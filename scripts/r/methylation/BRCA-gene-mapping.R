library(readr)
library(dplyr)
library(IlluminaHumanMethylation27k.db)

ProbeToSymbol <- IlluminaHumanMethylation27kSYMBOL
mapped_probes <- mappedkeys(ProbeToSymbol)
mapped_probes_List <- as.list(ProbeToSymbol[mapped_probes])

brca <- read_delim('../../KRAB_and_TCGA_Subtypes/Data v3 2016.06.06/results/BRCA/methylation.Clusters.csv', delim=',')

brca$probe <- unname(sapply(brca$gene, function(gene) {
  ns <- names(mapped_probes_List)[mapped_probes_List == gene]
  paste(unname(ns), collapse=' ')
}))

write.table(brca, file = 'methylation/results/BRCA-probes-mapped.csv', row.names = F, sep = ',')
